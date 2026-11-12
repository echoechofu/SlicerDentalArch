"""Mesial/distal width candidates for reviewed upper-arch crown components.

This is a semi-automatic geometric prototype, not a validated anatomical detector.
Its operational measurement is the largest surface-to-surface chord parallel to
an estimated tooth-local MD direction, within a reviewed crown search band.
Coordinates and widths are native STL units. Patient laterality is unverified.
"""
from pathlib import Path
import argparse, csv, hashlib, json
import numpy as np
from scipy.spatial import ConvexHull
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

from case_io import load_case, PREMOLARS

def read_stl(path,config_path):
    mesh=load_case(path,config_path)
    return tuple(mesh[k] for k in ('v','f','q','center','axes','face_groups'))

def polydata(v,f):
    points=vtk.vtkPoints();points.SetData(numpy_to_vtk(v,deep=True))
    cells=vtk.vtkCellArray();cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0,3*len(f)+1,3,dtype=np.int64),deep=True),numpy_to_vtkIdTypeArray(f.astype(np.int64).ravel(),deep=True))
    mesh=vtk.vtkPolyData();mesh.SetPoints(points);mesh.SetPolys(cells);return mesh

def unit(x):
    x=np.asarray(x,float);return x/np.linalg.norm(x)

def rotate(u,deg):
    a=np.deg2rad(deg);return np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])@u

def angle(u,v):return float(np.degrees(np.arctan2(u[0]*v[1]-u[1]*v[0],u@v)))

def manifold_face_components(f,face_ids):
    ff=f[face_ids]
    edges=np.concatenate([ff[:,[0,1]],ff[:,[1,2]],ff[:,[2,0]]]);owner=np.tile(np.arange(len(ff)),3)
    _,inverse,count=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    order=np.argsort(inverse,kind='stable');starts=np.r_[0,np.cumsum(count)][:-1][count==2]
    first=owner[order[starts]];second=owner[order[starts+1]]
    graph=coo_matrix((np.ones(len(first)),(first,second)),shape=(len(ff),len(ff))).tocsr()
    _,groups=connected_components(graph,directed=False)
    return groups

def tooth_region(q,f,labels,cfg):
    face_ids=np.flatnonzero(labels==cfg['group'])
    return np.unique(f[face_ids]),face_ids

def polygon_covariance(xy):
    # Area-uniform moments of the convex polygon avoid vertex-density bias.
    p=xy[ConvexHull(xy).vertices];r=np.roll(p,-1,axis=0)
    cross=p[:,0]*r[:,1]-r[:,0]*p[:,1];area=cross.sum()/2
    mean=((p+r)*cross[:,None]).sum(0)/(6*area)
    xx=(cross*(p[:,0]**2+p[:,0]*r[:,0]+r[:,0]**2)).sum()/(12*area)
    yy=(cross*(p[:,1]**2+p[:,1]*r[:,1]+r[:,1]**2)).sum()/(12*area)
    xy2=(cross*(2*p[:,0]*p[:,1]+p[:,0]*r[:,1]+r[:,0]*p[:,1]+2*r[:,0]*r[:,1])).sum()/(24*area)
    cov=np.array([[xx,xy2],[xy2,yy]])-np.outer(mean,mean)
    return cov,p

def estimate_axis(tq,cfg,prior):
    high=np.percentile(tq[:,2],99.5)
    if 'incisor' in cfg['kind']:
        upper=tq[tq[:,2]>np.percentile(tq[:,2],82),:2]
        cov,_=polygon_covariance(upper)
        vals,ev=np.linalg.eigh(cov);ix=int(np.argmax(np.abs(ev.T@prior)));u=ev[:,ix]
        method='Area-uniform PCA of upper 18% crown outline; tooth-order tangent selects/signs axis'
        metadata={'outline_eigenvalue_ratio':float(vals[1]/vals[0])}
    else:
        depth=4.0 if cfg['kind']=='molar' else 3.5
        band=tq[tq[:,2]>high-depth,:2]
        covariance,hull=polygon_covariance(band)
        # Orientation is selected from this tooth's outline, with a weak tooth-order
        # prior to disambiguate almost-square crowns. Do not maximize width over angles.
        trial=[]
        area0=np.ptp(hull@prior)*np.ptp(hull@rotate(prior,90))
        for deg in np.arange(-40,40.001,.5):
            u0=rotate(prior,deg);b0=rotate(u0,90)
            area=np.ptp(hull@u0)*np.ptp(hull@b0)
            score=area+.04*area0*(deg/30)**2
            trial.append((score,float(deg),float(area)))
        best=min(trial);u=rotate(prior,best[1]);method='Minimum-area upper-crown bounding rectangle with weak tooth-order orientation prior'
        near=[x[1] for x in trial if x[0]<=best[0]+.005*area0]
        vals=np.linalg.eigvalsh(covariance)
        metadata={'rectangle_angle_from_prior_deg':best[1],'rectangle_area':best[2],'orientation_search_limit_deg':40.,
                  'outline_eigenvalue_ratio':float(vals[1]/vals[0]),'axis_outline_depth':depth,
                  'near_optimal_angle_span_deg':float(max(near)-min(near))}
    if u@prior<0:u=-u
    return unit(u),method,metadata

def occlusal_height_map(local,faces,step):
    origin=local[:,:2].min(0)-step
    width,height=np.ceil((np.ptp(local[:,:2],axis=0)+2*step)/step).astype(int)+1
    h=np.full((height,width),-np.inf);face=np.full((height,width),-1,dtype=np.int64)
    for i,tri in enumerate(local[faces]):
        uv=(tri[:,:2]-origin)/step;x0,y0=np.maximum(np.floor(uv.min(0)).astype(int),0)
        x1,y1=np.minimum(np.ceil(uv.max(0)).astype(int),[width-1,height-1]);a,b,c=uv
        den=(b[1]-c[1])*(a[0]-c[0])+(c[0]-b[0])*(a[1]-c[1])
        if abs(den)<1e-10:continue
        xs,ys=np.meshgrid(np.arange(x0,x1+1),np.arange(y0,y1+1))
        p=((b[1]-c[1])*(xs-c[0])+(c[0]-b[0])*(ys-c[1]))/den
        r=((c[1]-a[1])*(xs-c[0])+(a[0]-c[0])*(ys-c[1]))/den
        z=p*tri[0,2]+r*tri[1,2]+(1-p-r)*tri[2,2]
        keep=(p>=-1e-8)&(r>=-1e-8)&(p+r<=1+1e-8)&(z>h[ys,xs])
        h[ys[keep],xs[keep]]=z[keep];face[ys[keep],xs[keep]]=i
    return h,face,origin

def central_groove_axis(q,f,ids,face_ids,cfg,prior,sigma=.18,step=.06):
    anchors=np.array(cfg['groove_anchors']);center=anchors.mean(0)
    b=unit(anchors[1]-anchors[0]);u=rotate(b,90)
    if u@prior<0:u=-u
    b=rotate(u,90);x=q[ids]
    loc=np.column_stack([(x[:,:2]-center)@u,(x[:,:2]-center)@b,x[:,2]])
    remap=np.full(len(q),-1);remap[ids]=np.arange(len(ids));tf=remap[f[face_ids]]
    h,face,origin=occlusal_height_map(loc,tf,step);valid=face>=0
    nearest=ndi.distance_transform_edt(~valid,return_distances=False,return_indices=True);filled=h[tuple(nearest)]
    sm=ndi.gaussian_filter(filled,sigma/step)
    curv=ndi.gaussian_filter(filled,sigma/step,order=(2,0))/step**2
    slope_v=ndi.gaussian_filter(filled,sigma/step,order=(1,0))/step
    ny,nx=h.shape;xs=origin[0]+np.arange(nx)*step;ys=origin[1]+np.arange(ny)*step
    mask=valid&(ndi.distance_transform_edt(valid)*step>.7)
    # A shape-initialized central corridor avoids gingival/side-wall minima.
    extent=np.linalg.norm(anchors[1]-anchors[0]);halfwidth=min(1.2,.4*extent)
    ulo,uhi=np.quantile(loc[:,0],[.25,.75]);cols=np.flatnonzero((xs>=ulo)&(xs<=uhi))
    normal_curvature=curv/(1+slope_v*slope_v)**1.5
    cost=-np.clip(normal_curvature,0,3)+.25*np.abs(slope_v)+.025*ys[:,None]**2
    allowed=mask&(abs(ys[:,None])<=halfwidth)
    cost=np.where(allowed,cost,np.inf)
    dp=cost[:,cols[0]]*step;back=[];delta=4
    for col in cols[1:]:
        stack=[]
        for shift in range(-delta,delta+1):
            tmp=np.roll(dp,shift)+.25*(shift*step)**2/step
            if shift>0:tmp[:shift]=np.inf
            elif shift<0:tmp[shift:]=np.inf
            stack.append(tmp)
        stack=np.array(stack);choice=stack.argmin(0);dp=stack[choice,np.arange(ny)]+cost[:,col]*step
        back.append(np.arange(ny)-(choice-delta))
    if not np.isfinite(dp).any():raise RuntimeError('Central groove corridor has no supported path; review anchors.')
    row=int(dp.argmin());rows=[row]
    for br in reversed(back):row=int(br[row]);rows.append(row)
    rows=np.array(rows[::-1]);uv=np.column_stack([xs[cols],ys[rows]])
    slope,intercept=np.polyfit(uv[:,0],uv[:,1],1);direction=unit(u+slope*b)
    points=[];barys=[];tri_ids=[]
    for row,col in zip(rows,cols):
        index=face[row,col];tri=loc[tf[index]]
        w=np.linalg.solve(np.vstack([tri[:,:2].T,np.ones(3)]),np.r_[xs[col],ys[row],1])
        points.append(w@q[f[face_ids[index]]]);barys.append(w.tolist());tri_ids.append(int(face_ids[index]))
    rms=float(np.std(uv[:,1]-slope*uv[:,0]-intercept));boundary=float(halfwidth-np.abs(uv[:,1]).max())
    evidence=np.maximum(normal_curvature[rows,cols],0)
    meta=dict(method='Central occlusal groove candidate: shape-initialized corridor, transverse concavity/slope, dynamic-programming continuity, fitted direction',
              anchors_aligned_xy=anchors.tolist(),anchors_are_provisional=True,corridor_half_width=halfwidth,
              anchor_perpendicular_axis_xy=u.tolist(),fitted_axis_xy=direction.tolist(),
              angle_from_anchor_axis_deg=angle(u,direction),line_fit_RMS=rms,
              corridor_boundary_distance=boundary,mean_transverse_concavity=float(evidence.mean()),
              median_transverse_concavity=float(np.median(evidence)),grid_step=step,smoothing_sigma=sigma,
              points_aligned=np.array(points).tolist(),triangle_indices=tri_ids,barycentric=barys,
              anatomical_groove_unverified=True,additional_identity_uncertainty=False,review_required=True)
    return direction,meta,dict(height=h,valid=valid,smoothed=sm,normal_curvature=normal_curvature,path_rows=rows,path_columns=cols,origin=origin,step=step)

def envelopes(local,faces,origin,step,shape):
    """Rasterize exact triangle intersections along MD onto a (BL,height) grid."""
    height,width=shape;lo=np.full(shape,np.inf);hi=np.full(shape,-np.inf)
    flo=np.full(shape,-1,np.int64);fhi=np.full(shape,-1,np.int64)
    for i,tri in enumerate(local[faces]):
        uv=(tri[:,1:]-origin)/step
        x0,y0=np.maximum(np.floor(uv.min(0)).astype(int),0)
        x1,y1=np.minimum(np.ceil(uv.max(0)).astype(int),[width-1,height-1])
        if x1<x0 or y1<y0:continue
        a,b,c=uv;den=(b[1]-c[1])*(a[0]-c[0])+(c[0]-b[0])*(a[1]-c[1])
        if abs(den)<1e-12:continue
        xs,ys=np.meshgrid(np.arange(x0,x1+1),np.arange(y0,y1+1))
        p=((b[1]-c[1])*(xs-c[0])+(c[0]-b[0])*(ys-c[1]))/den
        r=((c[1]-a[1])*(xs-c[0])+(a[0]-c[0])*(ys-c[1]))/den
        u=p*tri[0,0]+r*tri[1,0]+(1-p-r)*tri[2,0]
        valid=(p>=-1e-9)&(r>=-1e-9)&(p+r<=1+1e-9)
        less=valid&(u<lo[ys,xs]);more=valid&(u>hi[ys,xs])
        lo[ys[less],xs[less]]=u[less];flo[ys[less],xs[less]]=i
        hi[ys[more],xs[more]]=u[more];fhi[ys[more],xs[more]]=i
    return lo,hi,flo,fhi

def clipped_support(local,faces,face_ids,normals,bounds):
    """Exact extrema on eligible triangle patches clipped to the search rectangle.

    Also return the two surface locations; their Euclidean separation is NOT the
    supporting-plane width when their BL or height coordinates differ.
    """
    best=[-np.inf,np.inf];points=[None,None];tri_ids=[None,None];barys=[None,None]
    blo,bhi,zmin,zmax=bounds
    for j,face in enumerate(faces):
        side=0 if normals[j,0]>=.70 else (1 if normals[j,0]<=-.70 else -1)
        if side<0:continue
        tri=local[face]
        if tri[:,1].max()<blo or tri[:,1].min()>bhi or tri[:,2].max()<zmin or tri[:,2].min()>zmax:continue
        polygon=[(p,w) for p,w in zip(tri,np.eye(3))]
        for axis,threshold,sign in [(1,blo,1),(1,bhi,-1),(2,zmin,1),(2,zmax,-1)]:
            clipped=[]
            if not polygon:break
            previous,wp=polygon[-1];dp=sign*(previous[axis]-threshold)
            for current,wc in polygon:
                dc=sign*(current[axis]-threshold)
                if (dc>=0)!=(dp>=0):
                    t=dp/(dp-dc);clipped.append((previous+t*(current-previous),wp+t*(wc-wp)))
                if dc>=0:clipped.append((current,wc))
                previous,wp,dp=current,wc,dc
            polygon=clipped
        for p,w in polygon:
            better=p[0]>best[0] if side==0 else p[0]<best[1]
            if better:best[side]=float(p[0]);points[side]=p;tri_ids[side]=int(face_ids[j]);barys[side]=w
    if any(x is None for x in points):raise RuntimeError('No eligible proximal supporting patches.')
    return dict(width=best[0]-best[1],points_local=np.array(points),triangle_indices=tri_ids,barycentric=np.array(barys))

def find_pair(q,f,ids,face_ids,cfg,u,step=.10,depth_offset=0.,refine=True,component_face_ids=None):
    initial_step=step
    tq=q[ids];center=tq.mean(0);b=rotate(u,90)
    # b orientation does not affect maximization or MD labels.
    local=np.column_stack(((q[:,:2]-center[:2])@u,(q[:,:2]-center[:2])@b,q[:,2]))
    own=local[ids];zmax=float(np.percentile(own[:,2],99.5))
    depth=6.0 if ('incisor' in cfg['kind'] or cfg['kind']=='canine') else (5.0 if cfg['kind']=='premolar' else 4.5)
    exact_bottom=float(own[:,2].min());exact_top=float(own[:,2].max())
    posterior=cfg['kind'] in ['premolar','molar']
    hard_floor=exact_bottom+(exact_top-exact_bottom)*2/3
    zmin=hard_floor+max(depth_offset,0.) if posterior else zmax-depth+depth_offset
    band=own[own[:,2]>=zmin]
    blo,bhi=np.percentile(band[:,1],[10,90])
    zupper=zmax-.10
    origin=np.array([blo,zmin]);shape=(int(np.floor((zupper-zmin)/step))+1,int(np.floor((bhi-blo)/step))+1)
    ff=f[face_ids]
    triangles=local[ff]
    normals=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
    normals/=np.linalg.norm(normals,axis=1)[:,None]
    # Genuine opposing proximal faces must face the virtual caliper jaws.
    # This rejects occlusal and cervical cap surfaces in favor of proximal faces.
    proximal_mesial=normals[:,0]>=.70
    proximal_distal=normals[:,0]<=-.70
    def width_scores(lo,hi,flo,fhi):
        valid=(flo>=0)&(fhi>=0)
        valid &= proximal_mesial[np.maximum(fhi,0)]&proximal_distal[np.maximum(flo,0)]
        return np.where(valid,hi-lo,-np.inf)
    # Reject disjoint/intervening surfaces along the proposed chord using all
    # faces of the isolated tooth. Base surfaces are excluded by face adjacency.
    all_tri=local[f[component_face_ids if component_face_ids is not None else face_ids]]
    aa,bb,cc=all_tri[:,0,1:],all_tri[:,1,1:],all_tri[:,2,1:]
    denominator=(bb[:,1]-cc[:,1])*(aa[:,0]-cc[:,0])+(cc[:,0]-bb[:,0])*(aa[:,1]-cc[:,1])
    nonparallel=np.abs(denominator)>1e-12
    def crossings(yz):
        p=np.zeros(len(denominator));rr=np.zeros_like(p)
        p[nonparallel]=((bb[:,1]-cc[:,1])*(yz[0]-cc[:,0])+(cc[:,0]-bb[:,0])*(yz[1]-cc[:,1]))[nonparallel]/denominator[nonparallel]
        rr[nonparallel]=((cc[:,1]-aa[:,1])*(yz[0]-cc[:,0])+(aa[:,0]-cc[:,0])*(yz[1]-cc[:,1]))[nonparallel]/denominator[nonparallel]
        valid=nonparallel&(p>=-1e-8)&(rr>=-1e-8)&(p+rr<=1+1e-8)
        values=np.sort((p*all_tri[:,0,0]+rr*all_tri[:,1,0]+(1-p-rr)*all_tri[:,2,0])[valid])
        return values[np.r_[True,np.diff(values)>1e-5]] if len(values) else values
    rejections=0
    def select_valid(score,lo,hi,origin,spacing):
        nonlocal rejections
        order=np.argsort(score.ravel())[::-1]
        for ix in order:
            row,col=np.unravel_index(ix,score.shape)
            if not np.isfinite(score[row,col]):break
            yz=origin+np.array([col,row])*spacing
            hits=crossings(yz)
            extra=(hits>lo[row,col]+1e-5)&(hits<hi[row,col]-1e-5)
            if not extra.any():return row,col
            rejections+=1
        return None
    lo,hi,flo,fhi=envelopes(local,ff,origin,step,shape)
    score=width_scores(lo,hi,flo,fhi)
    if not np.isfinite(score).any():raise RuntimeError('No two-sided tooth intersections.')
    selected=select_valid(score,lo,hi,origin,step)
    if selected is None:raise RuntimeError('No eligible chord without additional intervening surfaces; tooth segmentation needs review.')
    row,col=selected;best=float(score[row,col]);yz=origin+np.array([col,row])*step
    ia,ib=int(fhi[row,col]),int(flo[row,col]);spacing=step
    if refine:
        for spacing in [step/5,step/25]:
            start=np.maximum(yz-step, [blo,zmin]);end=np.minimum(yz+step,[bhi,zupper])
            shape2=(int(np.floor((end[1]-start[1])/spacing))+1,int(np.floor((end[0]-start[0])/spacing))+1)
            l,h,fl,fh=envelopes(local,ff,start,spacing,shape2);s=width_scores(l,h,fl,fh)
            selected=select_valid(s,l,h,start,spacing)
            if selected is not None:
                rr2,cc2=selected
                if s[rr2,cc2]>best:
                    best=float(s[rr2,cc2]);yz=start+np.array([cc2,rr2])*spacing;ia,ib=int(fh[rr2,cc2]),int(fl[rr2,cc2])
            step=spacing
    pts=[];barys=[];native_faces=[]
    for ix in [ia,ib]:
        ti=face_ids[ix];tri=local[f[ti]]
        weights=np.linalg.solve(np.vstack([tri[:,1:].T,np.ones(3)]),np.r_[yz,1])
        pts.append(weights@q[f[ti]]);barys.append(weights);native_faces.append(int(ti))
    pts=np.array(pts)
    # Supporting-plane projection on the same clipped proximal patches is deliberately
    # separate from the chord and its Euclidean line measurement.
    support=clipped_support(local,ff,face_ids,normals,[blo,bhi,zmin,zupper])
    delta=pts[0]-pts[1]
    diag=dict(search_height_range=[zmin,zupper],search_BL_range=[float(blo),float(bhi)],
              selected_BL_height=yz.tolist(),boundary_height_distance=float(min(yz[1]-zmin,zupper-yz[1])),
              boundary_BL_distance=float(min(yz[0]-blo,bhi-yz[0])),
              proximal_patch_support_width=float(support['width']),support_minus_chord=float(support['width']-best),
              selected_face_normal_MD_components=[float(normals[ia,0]),float(normals[ib,0])],
              selected_face_normal_height_components=[float(normals[ia,2]),float(normals[ib,2])],
              rejected_rays_with_intervening_surfaces=int(rejections),
              full_component_intersections_along_selected_ray=crossings(yz).tolist(),
              direction_residual=float(np.linalg.norm(delta-np.r_[u,0]*np.linalg.norm(delta))),
              geometric_crown_bottom=exact_bottom,geometric_crown_top=exact_top,
              normalized_height_M_D=((pts[:,2]-exact_bottom)/(exact_top-exact_bottom)).tolist(),
              posterior_upper_third_enforced=posterior,posterior_hard_height_floor=hard_floor if posterior else None,
              final_sampling_step=spacing)
    support_q=np.array([w@q[f[ti]] for w,ti in zip(support['barycentric'],support['triangle_indices'])])
    support['points_aligned']=support_q;support.pop('points_local')
    return dict(points_aligned=pts,triangle_indices=native_faces,barycentric=np.array(barys),width=float(np.linalg.norm(delta)),diagnostics=diag,support=support),dict(score=score,origin=origin,initial_step=initial_step)

def run(source,out,step=.10,sensitivity=True,config_path=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if not .04<=step<=.20:raise ValueError('step must be between 0.04 and 0.20 native units.')
    mesh=load_case(source,config_path)
    v,f,q,center,axes,labels=[mesh[k] for k in ('v','f','q','center','axes','face_groups')]
    config=mesh['config'];region_configs=[dict(r) for r in config['regions']]
    anchors=config.get('groove_anchors',{})
    if any(name not in anchors for name in PREMOLARS):raise ValueError('All four premolar groove anchors are required; use run_all to derive provisional anchors from cusp candidates, or provide reviewed anchors.')
    for cfg in region_configs:
        if cfg['kind']=='premolar':cfg['groove_anchors']=anchors[cfg['id']]
    regions=[tooth_region(q,f,labels,cfg) for cfg in region_configs]
    centers=[]
    for ids,_ in regions:
        x=q[ids];centers.append(x[x[:,2]>np.percentile(x[:,2],99.5)-4,:2].mean(0))
    teeth=[]
    for i,cfg in enumerate(region_configs):
        ids,face_ids=regions[i]
        options={'component_face_ids':face_ids}
        if i%7==0:prior=centers[7 if i==0 else 0]-centers[i+1]
        elif i%7==6:prior=centers[i-1]-centers[i]
        else:prior=centers[i-1]-centers[i+1]
        prior=unit(prior);u,method,axis_meta=estimate_axis(q[ids],cfg,prior)
        if cfg['kind']=='premolar':
            outline_u=u.copy()
            u,groove,gdebug=central_groove_axis(q,f,ids,face_ids,cfg,prior)
            groove['points_native']=(np.array(groove['points_aligned'])@axes.T+center).tolist()
            np.savez_compressed(out/(cfg['id']+'_groove_search.npz'),**gdebug)
            variants=[]
            for sigma in [.14,.24]:
                direction,_,_=central_groove_axis(q,f,ids,face_ids,cfg,prior,sigma=sigma)
                variants.append(dict(sigma=sigma,direction_xy=direction.tolist(),angle_difference_deg=angle(u,direction)))
            groove['smoothing_sensitivity']=variants
            groove['max_smoothing_direction_change_deg']=max(abs(x['angle_difference_deg']) for x in variants)
            axis_meta={'groove':groove,'previous_outline_axis_xy':outline_u.tolist(),'angle_from_previous_outline_deg':angle(outline_u,u)}
            method=groove['method']
        pair,debug=find_pair(q,f,ids,face_ids,cfg,u,step,**options)
        np.savez_compressed(out/(cfg['id']+'_search.npz'),**debug)
        pts=pair.pop('points_aligned');barys=pair.pop('barycentric')
        support=pair['support'];support['points_native']=(support.pop('points_aligned')@axes.T+center).tolist();support['barycentric']=support['barycentric'].tolist()
        support['point_distance_3d']=float(np.linalg.norm(np.diff(support['points_native'],axis=0)))
        support['definition']='Projection width between two supporting planes on eligible proximal triangle patches; not the Euclidean distance of the support points.'
        native=pts@axes.T+center
        flags=[]
        if cfg['kind']=='premolar':
            g=axis_meta['groove']
            if g['additional_identity_uncertainty']:flags.append('central_groove_identity_requires_review')
            if g['max_smoothing_direction_change_deg']>5:flags.append('groove_direction_smoothing_sensitive')
            if g['corridor_boundary_distance']<.15:flags.append('groove_near_corridor_boundary')
            if g['line_fit_RMS']>.25:flags.append('groove_curved_or_branched')
        if abs(angle(prior,u))>25:flags.append('large_shape_axis_rotation_review')
        if axis_meta.get('near_optimal_angle_span_deg',0)>15:flags.append('outline_axis_ambiguous')
        d=pair['diagnostics']
        if d['boundary_height_distance']<.20:flags.append('near_height_search_boundary')
        if d['boundary_BL_distance']<.20:flags.append('near_BL_search_boundary')
        if d['support_minus_chord']>.5:flags.append('support_width_exceeds_chord_by_more_than_0.5')
        if d['rejected_rays_with_intervening_surfaces']>0:flags.append('overlapping_or_extra_surfaces_detected_in_search')
        tests={}
        if sensitivity:
            for deg in [-5,5]:
                p,_=find_pair(q,f,ids,face_ids,cfg,rotate(u,deg),step,refine=False,**options)
                tests['angle_%+d_deg'%deg]={'width':p['width'],'point_displacements':np.linalg.norm(p['points_aligned']-pts,axis=1).tolist()}
            for offset in ([.25,.5] if cfg['kind'] in ['premolar','molar'] else [-.5,.5]):
                p,_=find_pair(q,f,ids,face_ids,cfg,u,step,depth_offset=offset,refine=False,**options)
                tests['lower_height_%+.2f'%offset]={'actual_offset_native_units':offset,'width':p['width'],'point_displacements':np.linalg.norm(p['points_aligned']-pts,axis=1).tolist()}
            p,_=find_pair(q,f,ids,face_ids,cfg,u,step/2,**options)
            tests['half_step']={'width':p['width'],'point_displacements':np.linalg.norm(p['points_aligned']-pts,axis=1).tolist()}
            wd=max(abs(x['width']-pair['width']) for x in tests.values());pd=max(max(x['point_displacements']) for x in tests.values())
            if wd>.3:flags.append('width_parameter_sensitivity_over_0.3')
            if pd>1.:flags.append('point_parameter_sensitivity_over_1.0')
            d['max_tested_width_change']=float(wd);d['max_tested_point_change']=float(pd)
        teeth.append(dict(**cfg,**pair,mesial_native=native[0].tolist(),distal_native=native[1].tolist(),
            mesial_aligned=pts[0].tolist(),distal_aligned=pts[1].tolist(),barycentric=barys.tolist(),
            MD_axis_aligned=np.r_[u,0].tolist(),MD_axis_native=(np.r_[u,0]@axes.T).tolist(),
            arch_prior_axis_xy=prior.tolist(),axis_method=method,axis_diagnostics=axis_meta,
            segmentation_method='Reviewed whole-crown manifold-edge face component',
            tooth_face_count=int(len(face_ids)),
            axis_angle_from_arch_prior_deg=angle(prior,u),sensitivity=tests,review_flags=flags,review_required=True))
        print(cfg['id'],round(pair['width'],4),flags,flush=True)
    result=dict(source_file=str(Path(source).resolve()),source_sha256=config['sha256'],
        algorithm_version='2.0',method='Premolar groove-parallel chord; all posterior endpoints restricted to geometric occlusal third; anterior method preserved',
        units='Native STL coordinate units; nominal millimetres used for parameter selection, scale unconfirmed',
        coordinate_system='Raw STL coordinates; neither anatomical RAS/LPS nor patient laterality established',
        naming='VL/VR refer to the occlusal preview sides, not patient left/right. I1/I2/C/P1/P2/M1/M2 are reviewed tooth types. M points toward the anterior dental midline along each half-arch; D points away.',
        measurement_definition='max over permitted BL and height: u_max(BL,height)-u_min(BL,height), with opposing proximal normals. Premolar direction is fixed parallel to fitted central groove candidate in occlusal projection. All premolar/molar endpoints satisfy (z-z_min)/(z_max-z_min)>=2/3 on the isolated crown mesh; sensitivity trials cannot relax this bound. Anterior calculation is unchanged.',
        limitations=['Tooth-to-component labels and the supplied alignment require human review; source STL is unchanged.',
            'Fixed-direction same-height/BL chord approximates intrinsic mesiodistal crown width; it is not necessarily a physical contact between neighbors.',
            'Reference plane is a reviewed PCA occlusal approximation, not a clinically established occlusal plane or individual tooth axis.',
            'Posterior geometric upper-third constraint comes from user protocol, with no independently annotated CEJ.',
            'Premolar central groove candidates need anatomical review; smoothing sensitivity is not localization accuracy.',
            'Surface membership and parameter perturbations verify geometry/numerics, not anatomical ground truth.'],
        pca_center=center.tolist(),pca_axes_columns=axes.tolist(),initial_sampling_step=step,
        tooth_count=len(teeth),landmark_count=2*len(teeth),review_required=True,teeth=teeth)
    (out/'mesiodistal_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    with (out/'mesiodistal_measurements.csv').open('w',newline='',encoding='utf-8-sig') as file:
        writer=csv.writer(file);writer.writerow(['tooth_id','kind','M_x_native','M_y_native','M_z_native','D_x_native','D_y_native','D_z_native','candidate_width_native_units','M_triangle','D_triangle','normalized_height_M','normalized_height_D','upper_third_enforced','premolar_direction_change_deg','review_flags'])
        for t in teeth:writer.writerow([t['id'],t['kind'],*t['mesial_native'],*t['distal_native'],t['width'],*t['triangle_indices'],*t['diagnostics']['normalized_height_M_D'],t['diagnostics']['posterior_upper_third_enforced'],t['axis_diagnostics'].get('angle_from_previous_outline_deg',''),';'.join(t['review_flags'])])
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True);parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--config',type=Path,required=True);parser.add_argument('--step',type=float,default=.10);parser.add_argument('--skip-sensitivity',action='store_true')
    args=parser.parse_args()
    if not .04<=args.step<=.20:raise ValueError('--step must be between 0.04 and 0.20 native units.')
    run(args.input,args.output_dir,args.step,not args.skip_sensitivity,config_path=args.config)
