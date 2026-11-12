"""Multiscale cusp / central-fossa candidates for reviewed crown components.

Uses reviewed independent crown face components. Source triangles and native
coordinates are preserved; results are geometric candidates for review.
"""
from pathlib import Path
import argparse, json, csv
import numpy as np
from scipy import ndimage as ndi
from case_io import load_case, tooth_region
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray


def polydata(v, f):
    points=vtk.vtkPoints();points.SetData(numpy_to_vtk(v,deep=True))
    cells=vtk.vtkCellArray();cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0,3*len(f)+1,3,dtype=np.int64),deep=True),numpy_to_vtkIdTypeArray(f.astype(np.int64).ravel(),deep=True))
    mesh=vtk.vtkPolyData();mesh.SetPoints(points);mesh.SetPolys(cells);return mesh


def render(q,f,path,points=(),camera=(0,0,1),up=(0,1,0),size=(1400,1100),title=''):
    mesh=polydata(q,f)
    normals=vtk.vtkPolyDataNormals();normals.SetInputData(mesh);normals.SplittingOff();normals.ConsistencyOn();normals.AutoOrientNormalsOn()
    mapper=vtk.vtkPolyDataMapper();mapper.SetInputConnection(normals.GetOutputPort())
    actor=vtk.vtkActor();actor.SetMapper(mapper);actor.GetProperty().SetColor(.78,.80,.83);actor.GetProperty().SetAmbient(.25)
    ren=vtk.vtkRenderer();ren.SetBackground(.055,.075,.105);ren.AddActor(actor)
    for p in points:
        xyz=np.array(p['aligned_xyz']);col=(1,.28,.12) if p['kind']=='cusp_candidate' else (.1,.95,1)
        sphere=vtk.vtkSphereSource();sphere.SetCenter(*xyz);sphere.SetRadius(.20);sphere.SetThetaResolution(16);sphere.SetPhiResolution(16)
        mp=vtk.vtkPolyDataMapper();mp.SetInputConnection(sphere.GetOutputPort());ac=vtk.vtkActor();ac.SetMapper(mp);ac.GetProperty().SetColor(*col);ac.GetProperty().SetAmbient(.8);ren.AddActor(ac)
        tx=vtk.vtkBillboardTextActor3D();tx.SetInput(p['label']);tx.SetPosition(*(xyz+[.25,.25,.25]));tx.GetTextProperty().SetFontSize(16);tx.GetTextProperty().SetColor(*col);ren.AddActor(tx)
    if title:
        text=vtk.vtkTextActor();text.SetInput(title);text.GetTextProperty().SetFontSize(20);text.SetPosition(18,20);ren.AddViewProp(text)
    cam=ren.GetActiveCamera();mid=(q.min(0)+q.max(0))/2
    cam.SetFocalPoint(*mid);cam.SetPosition(*(mid+np.array(camera)*120));cam.SetViewUp(*up);cam.ParallelProjectionOn();ren.ResetCamera();cam.Zoom(1.10)
    win=vtk.vtkRenderWindow();win.SetOffScreenRendering(1);win.SetSize(*size);win.AddRenderer(ren);win.SetMultiSamples(4);win.Render()
    shot=vtk.vtkWindowToImageFilter();shot.SetInput(win);shot.SetInputBufferTypeToRGB();shot.ReadFrontBufferOff();shot.Update()
    wr=vtk.vtkPNGWriter();wr.SetFileName(str(path));wr.SetInputConnection(shot.GetOutputPort());wr.Write();win.Finalize()


def height_map(q,f,face_ids,step):
    x=q[np.unique(f[face_ids])];origin=x[:,:2].min(0)-step
    width,height=np.ceil((np.ptp(x[:,:2],axis=0)+2*step)/step).astype(int)+1
    h=np.full((height,width),-np.inf);face=np.full((height,width),-1,dtype=int)
    for fi in face_ids:
        tri=q[f[fi]];xy=(tri[:,:2]-origin)/step
        x0,y0=np.maximum(np.floor(xy.min(0)).astype(int),0);x1,y1=np.minimum(np.ceil(xy.max(0)).astype(int),[width-1,height-1])
        a,b,c=xy;den=(b[1]-c[1])*(a[0]-c[0])+(c[0]-b[0])*(a[1]-c[1])
        if abs(den)<1e-10:continue
        xs,ys=np.meshgrid(np.arange(x0,x1+1),np.arange(y0,y1+1))
        u=((b[1]-c[1])*(xs-c[0])+(c[0]-b[0])*(ys-c[1]))/den
        w=((c[1]-a[1])*(xs-c[0])+(a[0]-c[0])*(ys-c[1]))/den
        z=u*tri[0,2]+w*tri[1,2]+(1-u-w)*tri[2,2]
        keep=(u>=-1e-9)&(w>=-1e-9)&(u+w<=1+1e-9)&(z>h[ys,xs])
        h[ys[keep],xs[keep]]=z[keep];face[ys[keep],xs[keep]]=fi
    rows,cols=np.indices(h.shape);xy=np.stack([origin[0]+cols*step,origin[1]+rows*step],axis=-1)
    return h,face,xy


def detect(mesh,cfg,step=.08):
    q,f,v=mesh['q'],mesh['f'],mesh['v'];ids,face_ids=tooth_region(mesh,cfg)
    h,face,xy=height_map(q,f,face_ids,step);mask=face>=0
    dist=ndi.distance_transform_edt(mask)*step;interior=mask&(dist>1.0)
    if interior.sum()<3:raise ValueError(cfg['id']+': crown is too small for the native-unit detector parameters; check scale and tooth mapping.')
    origin=xy[interior].mean(0);design=np.column_stack([xy[interior]-origin,np.ones(interior.sum())])
    plane=np.linalg.lstsq(design,h[interior],rcond=None)[0]
    residual=h-((xy-origin)@plane[:2]+plane[2])
    near=ndi.distance_transform_edt(~mask,return_distances=False,return_indices=True);filled=residual[tuple(near)]
    scales=[.22,.38,.60];smooth=[ndi.gaussian_filter(filled,s/step) for s in scales];s=smooth[1]
    window=int(round(1.4/step))|1
    peaksets=[np.argwhere((sm>=ndi.maximum_filter(sm,size=window)-1e-9)&(dist>.55)) for sm in smooth]
    outer=ndi.gaussian_filter(filled,1.05/step);candidates=[]
    for row,col in peaksets[1]:
        persistence=sum(bool(len(p)) and float(np.min(np.linalg.norm(p-[row,col],axis=1)))*step<.8 for p in peaksets)
        prominence=float(s[row,col]-outer[row,col])
        if prominence<.1:continue
        score=prominence*(.6+.4*persistence/3)
        candidates.append(dict(score=score,row=int(row),col=int(col),scale_support=int(persistence),prominence=prominence))
    candidates.sort(key=lambda x:x['score'],reverse=True);selected=[]
    maxn={'canine':1,'premolar':2,'molar':4}[cfg['kind']]
    for item in candidates:
        if all(np.hypot(item['row']-a['row'],item['col']-a['col'])*step>1.8 for a in selected):selected.append(item)
        if len(selected)>=maxn:break

    def snap(item,kind,label):
        row,col=item['row'],item['col'];fi=int(face[row,col]);tri=q[f[fi]]
        w=np.linalg.solve(np.vstack([tri[:,:2].T,np.ones(3)]),np.r_[xy[row,col],1])
        native=w@v[f[fi]];aligned=w@tri
        flags=[]
        if item.get('scale_support',3)<3:flags.append('not_persistent_at_all_three_scales')
        if item.get('prominence',1)<.2:flags.append('low_prominence')
        return dict(label=label,region=cfg['id'],kind=kind,native_xyz=native.tolist(),aligned_xyz=aligned.tolist(),triangle_index=fi,barycentric=w.tolist(),
                    **{k:val for k,val in item.items() if k not in ['row','col']},review_flags=flags,review_required=True)
    points=[snap(item,'cusp_candidate',cfg['id']+'_C'+str(k)) for k,item in enumerate(selected,1)]
    pool=[]
    if cfg['kind']=='molar':
        ctr=np.median(xy[mask],axis=0);fx=ndi.gaussian_filter(filled,.30/step)
        loc=(fx<=ndi.minimum_filter(fx,size=int(round(1.2/step))|1)+1e-9)&(dist>1.35)
        for row,col in np.argwhere(loc):
            central=float(np.linalg.norm(xy[row,col]-ctr))
            if central>3.:continue
            depth=float(outer[row,col]-fx[row,col])
            if depth<.06:continue
            hxx=(fx[row,col+2]-2*fx[row,col]+fx[row,col-2])/(2*step)**2
            hyy=(fx[row+2,col]-2*fx[row,col]+fx[row-2,col])/(2*step)**2
            hxy=(fx[row+2,col+2]-fx[row+2,col-2]-fx[row-2,col+2]+fx[row-2,col-2])/(16*step**2)
            eig=np.linalg.eigvalsh([[hxx,hxy],[hxy,hyy]])
            if eig[0]<-.02:continue
            score=depth*np.exp(-.5*(central/2.)**2)
            pool.append(dict(score=float(score),row=int(row),col=int(col),local_depth=depth,distance_to_roi_center=central,hessian_eigenvalues=eig.tolist()))
        pool.sort(key=lambda x:x['score'],reverse=True)
        if pool:points.append(snap(pool[0],'central_fossa_candidate',cfg['id']+'_F'))
    debug=dict(region=cfg['id'],kind=cfg['kind'],cusp_pool=[snap(p,'cusp_candidate',cfg['id']+'_P'+str(i)) for i,p in enumerate(candidates)],fossa_pool=[snap(p,'central_fossa_candidate',cfg['id']+'_FP'+str(i)) for i,p in enumerate(pool)],cusp_expected_cap=maxn,cusp_selected=len(selected))
    return points,debug


def run(out,source=None,config_path=None,step=.08,previews=True):
    if not .04<=step<=.16:raise ValueError('step must be between 0.04 and 0.16 native units.')
    out=Path(out);out.mkdir(parents=True,exist_ok=True);mesh=load_case(source,config_path)
    result=[];debug=[]
    for cfg in mesh['config']['regions']:
        if 'incisor' in cfg['kind']:continue
        points,d=detect(mesh,cfg,step);result.extend(points);debug.append(d)
        if previews:
            ids,ff=tooth_region(mesh,cfg);remap=np.full(len(mesh['q']),-1);remap[ids]=np.arange(len(ids))
            render(mesh['q'][ids],remap[mesh['f'][ff]],out/(cfg['id']+'_top.png'),points,size=(650,600),title=cfg['id'])
            render(mesh['q'][ids],remap[mesh['f'][ff]],out/(cfg['id']+'_oblique.png'),points,camera=(.25,-1,1),up=(0,0,1),size=(650,600),title=cfg['id'])
    data=dict(source_sha256=mesh['config']['sha256'],units='native STL coordinate units',algorithm='Per-crown detrended multiscale maxima / central interior concavity; adapted reviewed crown map',parameters=dict(grid_step=step,scales=[.22,.38,.60],prominence_min=.1,cusp_separation=1.8),points=result,
        limitations=['Geometric candidates, no manual gold standard; peak labels are score-ranked and not named anatomical cusps.','Upper second molars may have fewer discernible cusps; detector caps counts but does not force missing candidates.','VL/VR are preview sides, not patient laterality.'],review_required=True)
    (out/'landmark_candidates.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf8')
    (out/'candidates.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    (out/'candidate_debug.json').write_text(json.dumps(debug,indent=2),encoding='utf8')
    with (out/'landmark_candidates.csv').open('w',encoding='utf-8-sig',newline='') as file:
        w=csv.writer(file);w.writerow(['label','region','kind','native_x','native_y','native_z','score','triangle_index','review_flags'])
        for p in result:w.writerow([p['label'],p['region'],p['kind'],*p['native_xyz'],p['score'],p['triangle_index'],';'.join(p['review_flags'])])
    if previews:render(mesh['q'],mesh['f'],out/'landmark_preview.png',result,size=(1700,1400))
    print(json.dumps({'points':len(result),'by_tooth':{x['region']:x['cusp_selected'] for x in debug}},indent=2),flush=True)
    return data


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--input',type=Path,required=True);p.add_argument('--config',type=Path,required=True);p.add_argument('--step',type=float,default=.08);p.add_argument('--no-previews',action='store_true');a=p.parse_args()
    run(a.output_dir,a.input,a.config,a.step,not a.no_previews)
