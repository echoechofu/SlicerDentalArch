"""Incisal-edge candidates from four reviewed incisor crown components.
Uses reviewed tooth-face groups, local multiscale ridges, continuous paths,
and projection onto the source STL. No canine ridges are included.
"""
from pathlib import Path
import argparse, csv, json
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
from case_io import load_case, tooth_region
def polydata(v,f):
    points=vtk.vtkPoints();points.SetData(numpy_to_vtk(v,deep=True))
    cells=vtk.vtkCellArray();cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0,3*len(f)+1,3,dtype=np.int64),deep=True),numpy_to_vtkIdTypeArray(f.astype(np.int64).ravel(),deep=True))
    mesh=vtk.vtkPolyData();mesh.SetPoints(points);mesh.SetPolys(cells);return mesh

def height_map(local,faces,step):
    origin=local[:,:2].min(0)-step
    width,height=np.ceil((np.ptp(local[:,:2],axis=0)+2*step)/step).astype(int)+1
    h=np.full((height,width),-np.inf);face=np.full((height,width),-1,dtype=np.int64)
    for i,tri in enumerate(local[faces]):
        uv=(tri[:,:2]-origin)/step
        x0,y0=np.maximum(np.floor(uv.min(0)).astype(int),0)
        x1,y1=np.minimum(np.ceil(uv.max(0)).astype(int),[width-1,height-1])
        xs,ys=np.meshgrid(np.arange(x0,x1+1),np.arange(y0,y1+1));a,b,c=uv
        den=(b[1]-c[1])*(a[0]-c[0])+(c[0]-b[0])*(a[1]-c[1])
        if abs(den)<1e-10:continue
        p=((b[1]-c[1])*(xs-c[0])+(c[0]-b[0])*(ys-c[1]))/den
        r=((c[1]-a[1])*(xs-c[0])+(a[0]-c[0])*(ys-c[1]))/den
        z=p*tri[0,2]+r*tri[1,2]+(1-p-r)*tri[2,2]
        keep=(p>=-1e-8)&(r>=-1e-8)&(p+r<=1+1e-8)&(z>h[ys,xs])
        h[ys[keep],xs[keep]]=z[keep];face[ys[keep],xs[keep]]=i
    return h,face,origin

def closest_project(mesh,coordinates):
    locator=vtk.vtkStaticCellLocator();locator.SetDataSet(mesh);locator.BuildLocator()
    out=[];faces=[];distances=[]
    for p in coordinates:
        result=[0.,0.,0.];cell=vtk.reference(0);sub=vtk.reference(0);dist=vtk.reference(0.)
        locator.FindClosestPoint(p,result,cell,sub,dist)
        out.append(result);faces.append(int(cell));distances.append(float(dist)**.5)
    return np.array(out),np.array(faces),np.array(distances)

def extract_curve(v,f,q,ids,face_ids,step=.08,smoothing=.16,arch_center=None):
    if arch_center is None or np.asarray(arch_center).shape!=(2,) or not np.isfinite(arch_center).all():
        raise ValueError('A finite aligned XY arch-center prior is required.')
    remap=np.full(len(v),-1,dtype=int);remap[ids]=np.arange(len(ids));tf=remap[f[face_ids]];tq=q[ids]
    tooth_center=tq.mean(0)
    # PCA of the upper crown band estimates the transverse direction; the reviewed arch
    # tangent is used to choose an axis instead of confusing crown depth with edge width.
    band=tq[tq[:,2]>np.percentile(tq[:,2],82)]
    _,ab=np.linalg.eigh(np.cov(band[:,:2].T))
    radial=tooth_center[:2]-arch_center;radial/=np.linalg.norm(radial)
    expected=np.array([radial[1],-radial[0]])
    u=ab[:,np.argmax(np.abs(ab.T@expected))]
    if u@expected<0:u=-u
    labial=np.array([-u[1],u[0]])
    if labial@radial<0:labial=-labial
    local=np.column_stack([(tq[:,:2]-tooth_center[:2])@u,(tq[:,:2]-tooth_center[:2])@labial,tq[:,2]])
    h,face,origin=height_map(local,tf,step)
    valid=face>=0
    near=ndi.distance_transform_edt(~valid,return_distances=False,return_indices=True)
    filled=h[tuple(near)]
    sm=ndi.gaussian_filter(filled,smoothing/step)
    sm2=ndi.gaussian_filter(filled,.28/step)
    # Prefer a high, convex ridge on the labial half of the crown; suppress cingulum.
    rows,cols=np.indices(h.shape);uv=np.stack([origin[0]+cols*step,origin[1]+rows*step],axis=-1)
    zmax=np.max(np.where(valid,sm,-np.inf),axis=0)
    supported=np.any(valid,axis=0)
    finite_cols=np.flatnonzero(supported)
    middle=finite_cols[len(finite_cols)//4:3*len(finite_cols)//4]
    trend=np.polyfit(origin[0]+middle*step,zmax[middle],1)
    incisal_baseline=np.polyval(trend,origin[0]+np.arange(h.shape[1])*step)
    eligible=supported&(zmax>=incisal_baseline-.65)
    # Keep a contiguous extent around the center; discard rounded side walls.
    regions,num=ndi.label(eligible)
    center_col=int(np.median(middle));region=regions[center_col]
    if region==0:
        region=max(range(1,num+1),key=lambda r:np.sum(regions==r))
    columns=np.flatnonzero(regions==region)
    trim=max(1,int(round(.16/step)))
    columns=columns[trim:-trim]
    if len(columns)<10:raise RuntimeError('Insufficient supported incisal sections')
    slope=ndi.gaussian_filter(filled,smoothing/step,order=(1,0))/step
    convex=np.maximum(-ndi.gaussian_filter(filled,smoothing/step,order=(2,0))/step**2,0)/(1+slope**2)**1.5
    convex2=np.maximum(-ndi.gaussian_filter(filled,.28/step,order=(2,0))/step**2,0)
    # Height provides the main constraint; curvature only disambiguates near-top bands.
    cost=np.maximum(zmax[None,:]-sm,0)+.10*np.minimum(np.abs(sm-sm2),1)-.045*np.clip(convex,0,3)-.02*np.clip(convex2,0,3)
    allowed=valid&(sm>incisal_baseline[None,:]-1.0)&(uv[:,:,1]>np.percentile(local[:,1],45))
    cost=np.where(allowed,cost,np.inf)
    # Dynamic programming over transverse sections penalizes jumps in labiolingual position.
    delta=int(np.ceil(.40/step));penalty=.052
    # Discretize integral unary_cost du + lambda (dv/du)^2 du consistently.
    dp=cost[:,columns[0]].copy()*step;back=[]
    for col in columns[1:]:
        stacks=[]
        for shift in range(-delta,delta+1):
            candidate=np.roll(dp,shift)+penalty*(shift*step)**2/step
            if shift>0:candidate[:shift]=np.inf
            if shift<0:candidate[shift:]=np.inf
            stacks.append(candidate)
        stack=np.array(stacks);choice=stack.argmin(0)
        dp=cost[:,col]*step+stack[choice,np.arange(len(dp))]
        back.append(np.arange(len(dp))-(choice-delta))
    if not np.isfinite(dp).any():raise RuntimeError('No continuous supported path')
    ridx=int(np.argmin(dp));path=[ridx]
    for b in reversed(back):ridx=int(b[ridx]);path.append(ridx)
    path=np.array(path[::-1]);labial_path=origin[1]+path*step
    labial_path=ndi.gaussian_filter1d(labial_path,.18/step,mode='nearest')
    target=[]
    for col,lcoord in zip(columns,labial_path):
        row=int(np.clip(round((lcoord-origin[1])/step),0,h.shape[0]-1))
        if not valid[row,col]:row=int(path[len(target)])
        uv_target=np.array([origin[0]+col*step,origin[1]+row*step])
        fi=face[row,col];tri=local[tf[fi]]
        bary=np.linalg.solve(np.vstack([tri[:,:2].T,np.ones(3)]),np.r_[uv_target,1])
        target.append(bary@v[ids][tf[fi]])
    curve=np.array(target)
    # Smooth the path mildly, then reproject every control point to this tooth only.
    smooth_curve=ndi.gaussian_filter1d(curve,.12/step,axis=0,mode='nearest')
    mesh=polydata(v[ids],tf);curve,local_face,snap_distance=closest_project(mesh,smooth_curve)
    length=np.linalg.norm(np.diff(curve,axis=0),axis=1);arc=np.r_[0,np.cumsum(length)]
    mid=np.array([np.interp(arc[-1]/2,arc,curve[:,j]) for j in range(3)])
    mid,midface,_=closest_project(mesh,[mid]);mid=mid[0]
    # Surface distance at segment midpoints reports chord approximation, not anatomical accuracy.
    _,_,chord_error=closest_project(mesh,(curve[:-1]+curve[1:])/2)
    barycentric=[]
    for p,fi in zip(curve,local_face):
        tri=v[ids][tf[fi]]
        weights=np.linalg.lstsq(np.vstack([tri.T,np.ones(3)]),np.r_[p,1],rcond=None)[0]
        barycentric.append(weights.tolist())
    return dict(points_native=curve.tolist(),triangle_indices=face_ids[local_face].tolist(),barycentric=barycentric,
                endpoint_A_native=curve[0].tolist(),endpoint_B_native=curve[-1].tolist(),midpoint_native=mid.tolist(),
                arc_length_native_units=float(arc[-1]),local_transverse_axis_xy=u.tolist(),local_labial_axis_xy=labial.tolist(),
                max_segment_midpoint_surface_distance=float(chord_error.max()),max_projection_displacement=float(snap_distance.max()),
                control_point_count=len(curve),review_required=True),dict(local=local,height=h,smoothed=sm,face=face,columns=columns,path_rows=path,origin=origin)


def dense(curve,n=800):
    p=np.array(curve['points_native']);a=np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))]
    return np.column_stack([np.interp(np.linspace(0,a[-1],n),a,p[:,j]) for j in range(3)])


def run(source=None,out=None,config_path=None,step=.08,smoothing=.16,sensitivity=True):
    if out is None:raise ValueError('An output directory is required.')
    if not .04<=step<=.16 or not .08<=smoothing<=.32:raise ValueError('Sampling/smoothing outside the supported native-unit parameter range.')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    mesh=load_case(source,config_path);v,f,q,center,axes=[mesh[k] for k in ['v','f','q','center','axes']]
    config=mesh['config'];regions={r['id']:r for r in config['regions']}
    # The radial prior is derived from the two reviewed first-molar crowns.
    arch_center=np.mean([q[tooth_region(mesh,regions[k])[0],:2].mean(0) for k in ['VL_M1','VR_M1']],axis=0)
    curves=[];checks=[];sensitivity_rows=[]
    for name,region in [('E1','VL_I2'),('E2','VL_I1'),('E3','VR_I1'),('E4','VR_I2')]:
        ids,face_ids=tooth_region(mesh,regions[region])
        curve,debug=extract_curve(v,f,q,ids,face_ids,step,smoothing,arch_center)
        curve.update(id=name,region=region,face_group=regions[region]['group'])
        curve['review_flags']=['algorithmic_edge_extent_needs_anatomical_review']
        p=np.array(curve['points_native']);w=np.array(curve['barycentric']);fi=np.array(curve['triangle_indices'])
        reconstructed=np.einsum('ni,nij->nj',w,v[f[fi]])
        error=np.linalg.norm(p-reconstructed,axis=1)
        assert np.all(np.isin(fi,face_ids)) and error.max()<1e-6
        assert w.min()>-1e-6 and w.max()<1+1e-6 and np.max(abs(w.sum(1)-1))<1e-6
        surface=polydata(v,f[face_ids])
        _,_,ld=closest_project(surface,[curve[k] for k in ['endpoint_A_native','midpoint_native','endpoint_B_native']])
        samples=np.concatenate([(1-t)*p[:-1]+t*p[1:] for t in [.1,.25,.5,.75,.9]])
        _,_,sd=closest_project(surface,samples)
        assert ld.max()<1e-6
        checks.append(dict(id=name,region=region,control_points=len(p),max_barycentric_error=float(error.max()),min_barycentric_weight=float(w.min()),max_landmark_surface_distance=float(ld.max()),max_sampled_segment_surface_distance=float(sd.max()),all_points_on_assigned_crown=True))
        if sensitivity:
            for trial_step in [.06,.10]:
                trial,_=extract_curve(v,f,q,ids,face_ids,trial_step,smoothing,arch_center)
                a,b=dense(curve),dense(trial);dist=np.r_[cKDTree(a).query(b)[0],cKDTree(b).query(a)[0]]
                delta_mid=float(np.linalg.norm(np.array(curve['midpoint_native'])-trial['midpoint_native']))
                delta_end=max(float(np.linalg.norm(np.array(curve[k])-trial[k])) for k in ['endpoint_A_native','endpoint_B_native'])
                sensitivity_rows.append(dict(id=name,region=region,step=trial_step,midpoint_shift=delta_mid,max_endpoint_shift=delta_end,sampled_symmetric_curve_distance_max=float(dist.max()),arc_length_change=trial['arc_length_native_units']-curve['arc_length_native_units']))
                if max(delta_mid,delta_end)>.4:curve['review_flags'].append('sampling_sensitive_extent')
        curves.append(curve);np.savez_compressed(out/(name+'_debug.npz'),**debug)
    result=dict(source_file=config['source'],source_sha256=config['sha256'],method='Reviewed face groups + original multiscale incisal ridge/dynamic-programming algorithm + original-surface projection',
        units='Native STL coordinate units, no physical-unit assumption',coordinate_system='Unchanged source STL coordinates; VL/VR refer to preview, not patient laterality',scope='Four incisors only; no canine ridges',
        grid_spacing=step,smoothing_sigma=smoothing,pca_center=center.tolist(),pca_axes_columns=axes.tolist(),arch_center_prior_aligned_xy=arch_center.tolist(),
        labels='E1=VL_I2, E2=VL_I1, E3=VR_I1, E4=VR_I2. A/B are truncated curve endpoints, not assigned mesial/distal contacts.',
        midpoint_definition='Half polyline arc length followed by projection onto the same tooth surface',review_required=True,
        limitations=['Case-specific reviewed face-group labels; not a universal tooth classifier.','Endpoint extent retains the existing 0.65-unit height-drop criterion, so anatomical corners need review.','Surface projection and numerical stability do not establish anatomical accuracy.','Worn or bevelled incisal surfaces can admit more than one representative edge.'],curves=curves)
    (out/'incisal_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    (out/'validation.json').write_text(json.dumps(dict(geometry=checks,sampling_sensitivity=sensitivity_rows,source_unchanged=True,note='Numerical geometry and sensitivity only; not clinical localization accuracy.'),indent=2),encoding='utf8')
    with (out/'incisal_control_points.csv').open('w',newline='',encoding='utf-8-sig') as file:
        writer=csv.writer(file);writer.writerow(['curve_id','tooth_id','point_index','x_native','y_native','z_native','source_triangle_index'])
        for curve in curves:
            for i,(p,fi) in enumerate(zip(curve['points_native'],curve['triangle_indices'])):writer.writerow([curve['id'],curve['region'],i,*p,fi])
    with (out/'incisal_landmarks.csv').open('w',newline='',encoding='utf-8-sig') as file:
        writer=csv.writer(file);writer.writerow(['curve_id','tooth_id','landmark','x_native','y_native','z_native'])
        for curve in curves:
            for key in ['endpoint_A_native','midpoint_native','endpoint_B_native']:writer.writerow([curve['id'],curve['region'],key,*curve[key]])
    print(json.dumps([dict(id=c['id'],tooth=c['region'],points=c['control_point_count'],length=c['arc_length_native_units'],flags=c['review_flags']) for c in curves],indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True);parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True);parser.add_argument('--step',type=float,default=.08)
    parser.add_argument('--smoothing',type=float,default=.16);parser.add_argument('--skip-sensitivity',action='store_true')
    args=parser.parse_args()
    if not .04<=args.step<=.16 or not .08<=args.smoothing<=.32:raise ValueError('Sampling/smoothing outside reviewed parameter range')
    run(args.input,args.output_dir,args.config,args.step,args.smoothing,not args.skip_sensitivity)


