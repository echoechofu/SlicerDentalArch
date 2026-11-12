"""Current-form reference arch from accepted V2 MD observations.

Planar, asymmetric, endpoint-constrained robust cubic B-spline approximation.
This is an explicit measurement-reference protocol, not an ideal treatment arch.
"""
from pathlib import Path
import json,csv,argparse,hashlib
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize
from scipy.integrate import quad
from scipy.spatial import cKDTree

TOOTH_ORDER=['VL_P2','VL_P1','VL_C','VL_I2','VL_I1','VR_I1','VR_I2','VR_C','VR_P1','VR_P2']

def observations(data):
    byid={t['id']:t for t in data['teeth']}
    if data.get('algorithm_version')!='2.0' or not data.get('source_sha256'):
        raise ValueError('Expected V2 measurements with a source SHA-256 digest.')
    required=set(TOOTH_ORDER+['VL_M1','VR_M1'])
    if not required.issubset(byid) or len(byid)!=len(data['teeth']):raise ValueError('Missing or duplicate tooth measurements for reference-arch fitting.')
    center=np.array(data['pca_center']);axes=np.array(data['pca_axes_columns'])
    if center.shape!=(3,) or axes.shape!=(3,3) or not np.isfinite(center).all() or not np.isfinite(axes).all():raise ValueError('Invalid source frame.')
    if not np.allclose(axes.T@axes,np.eye(3),atol=1e-10,rtol=0) or not np.isclose(np.linalg.det(axes),1,atol=1e-10,rtol=0):raise ValueError('Non-orthonormal source frame.')
    out=[]
    for name in TOOTH_ORDER:
        for key in (['distal','mesial'] if name.startswith('VL') else ['mesial','distal']):
            t=byid[name];point=np.array(t[key+'_native']);q=(point-center)@axes
            out.append(dict(id=name+'_'+key[0].upper(),tooth=name,role=key,native=point.tolist(),aligned=q.tolist(),weight=.5,
                source_triangle=t['triangle_indices'][0 if key=='mesial' else 1],source_review_flags=t['review_flags']))
    endpoints=np.array([byid['VL_M1']['mesial_native'],byid['VR_M1']['mesial_native']]);eq=(endpoints-center)@axes
    if not np.isfinite(endpoints).all() or not all(np.isfinite(o['native']).all() for o in out):raise ValueError('Non-finite measurement coordinates.')
    return out,endpoints,eq,center,axes

def basis(n=7):
    degree=3;interior=np.linspace(0,1,n-degree+1)[1:-1]
    knots=np.r_[np.zeros(4),interior,np.ones(4)]
    return BSpline(knots,np.eye(n),degree),knots

def bending_matrix(bs):
    k=np.zeros((bs.c.shape[0],bs.c.shape[0]));nodes,weights=np.polynomial.legendre.leggauss(3)
    breaks=np.unique(bs.t)
    for a,b in zip(breaks[:-1],breaks[1:]):
        t=(nodes+1)*(b-a)/2+a;B=bs(t,nu=2);k+=B.T@(B*(weights*(b-a)/2)[:,None])
    return k

def footpoints(curve,points):
    # Enumerate all stationary points of squared Euclidean distance within each
    # cubic span. The derivative is degree 5, so global nearest footpoints do not
    # depend on a single local optimizer or an initial correspondence.
    breaks=np.unique(curve.t);ts=[];dist=[]
    for p in points:
        candidates=list(breaks)
        for a,b in zip(breaks[:-1],breaks[1:]):
            # Use normalized span coordinate s=(t-a)/(b-a). BSpline evaluates
            # right-hand span derivatives at interior knots. Normalization and
            # trimming numerical leading zeroes protect linear/quadratic spans.
            span=b-a
            coeff=np.array([curve(a),curve(a,nu=1)*span,curve(a,nu=2)*span**2/2,curve(a,nu=3)*span**3/6])
            pol=np.zeros(6)
            for axis in range(2):
                co=coeff[:,axis].copy();co[0]-=p[axis]
                piece=np.polynomial.polynomial.polymul(co,np.array([coeff[1,axis],2*coeff[2,axis],3*coeff[3,axis]]))
                pol[:len(piece)]+=piece
            cutoff=1e-12*np.max(np.abs(pol))
            while len(pol)>1 and abs(pol[-1])<=cutoff:pol=pol[:-1]
            roots=np.polynomial.polynomial.polyroots(pol)
            for root in roots:
                if abs(root.imag)<1e-7:
                    t=float(a+root.real*span)
                    if a<t<b:candidates.append(t)
        candidates=np.array(candidates);d2=np.sum((curve(candidates)-p)**2,axis=1);i=int(np.argmin(d2))
        ts.append(candidates[i]);dist.append(float(np.sqrt(d2[i])))
    return np.array(ts),np.array(dist)

def huber_values(r,delta):return np.where(r<=delta,.5*r*r,delta*(r-.5*delta))

def fit_arch(points,endpoints,weights=None,n_control=7,smooth=1e-4,huber_fraction=.015,max_iter=300):
    points=np.asarray(points,float);ends=np.asarray(endpoints,float)
    weights=np.full(len(points),.5) if weights is None else np.asarray(weights,float)
    origin=ends.mean(0);scale=float(np.linalg.norm(ends[1]-ends[0]));p=(points-origin)/scale;e=(ends-origin)/scale
    bs,knots=basis(n_control);K=bending_matrix(bs)
    sequence=np.vstack([e[0],p,e[1]]);steps=np.linalg.norm(np.diff(sequence,axis=0),axis=1)
    t=np.cumsum(steps)[:-1]/sum(steps);robust=np.ones(len(points));control=np.zeros((n_control,2));control[[0,-1]]=e
    free=np.arange(1,n_control-1);fixed=np.array([0,n_control-1]);history=[];last=None;converged=False
    for it in range(25):
        B=bs(t);w=weights*robust;normal=B.T@(w[:,None]*B)+smooth*K
        rhs=B.T@(w[:,None]*p)-normal[:,fixed]@e
        control[free]=np.linalg.solve(normal[np.ix_(free,free)],rhs[free])
        curve=BSpline(knots,control,3);new_t,dist=footpoints(curve,p)
        new_robust=np.minimum(1,huber_fraction/np.maximum(dist,1e-14))
        objective=float(np.dot(weights,huber_values(dist,huber_fraction))+.5*smooth*np.sum(control*(K@control)))
        sample=curve(np.linspace(0,1,101))
        change=float(np.max(np.linalg.norm(sample-last,axis=1))) if last is not None else np.inf
        history.append(dict(iteration=it+1,objective=objective,max_curve_change_native=change*scale if np.isfinite(change) else None))
        t=new_t;robust=new_robust;last=sample
        if change<1e-9 and it>3:break
    # Refine the controls directly. Closest-point minimization supplies footpoints;
    # the envelope theorem gives the exact first derivative with respect to C.
    # This avoids the very slow terminal convergence of alternating projection.
    def objective_gradient(flat):
        ctrl=control.copy();ctrl[free]=flat.reshape(-1,2);candidate=BSpline(knots,ctrl,3)
        ft,ds=footpoints(candidate,p);B=bs(ft);res=candidate(ft)-p
        rw=np.minimum(1,huber_fraction/np.maximum(ds,1e-14))*weights
        val=float(np.dot(weights,huber_values(ds,huber_fraction))+.5*smooth*np.sum(ctrl*(K@ctrl)))
        grad=B.T@(rw[:,None]*res)+smooth*(K@ctrl)
        return val,grad[free].ravel()
    optimized=minimize(objective_gradient,control[free].ravel(),jac=True,method='L-BFGS-B',
        options={'maxiter':max_iter,'ftol':1e-14,'gtol':1e-9,'maxls':40,'maxcor':15})
    control[free]=optimized.x.reshape(-1,2)
    objective,gradient=objective_gradient(optimized.x)
    converged=bool(optimized.success and np.linalg.norm(gradient,ord=np.inf)<1e-7)
    native_control=control*scale+origin;curve=BSpline(knots,native_control,3)
    t,dist=footpoints(curve,points);length,quad_error=quad(lambda x:float(np.linalg.norm(curve(x,nu=1))),0,1,points=np.unique(knots)[1:-1],epsabs=1e-10,epsrel=1e-10,limit=100)
    sample_t=np.linspace(0,1,1001);samples=curve(sample_t);v=curve(sample_t,nu=1);acc=curve(sample_t,nu=2)
    cross=v[:,0]*acc[:,1]-v[:,1]*acc[:,0];curvature=cross/np.linalg.norm(v,axis=1)**3
    return dict(curve=curve,control=native_control,knots=knots,t=t,residuals=dist,robust_weights=np.minimum(1,huber_fraction*scale/np.maximum(dist,1e-14)),
        length=float(length),integration_error=float(quad_error),history=history,converged=converged,iterations=int(optimized.nit),
        optimizer=dict(method='25 IRLS/projection warm-up steps then L-BFGS-B on exact closest-point Huber objective',success=bool(optimized.success),message=str(optimized.message),objective=objective,gradient_infinity_norm=float(np.linalg.norm(gradient,ord=np.inf))),
        width_scale=scale,huber_delta_native=huber_fraction*scale,shape_diagnostics=dict(min_speed=float(np.linalg.norm(v,axis=1).min()),
        min_signed_curvature=float(curvature.min()),max_signed_curvature=float(curvature.max()),
        transverse_extent=float(np.ptp(samples[:,0])),depth=float(samples[:,1].max()-ends[:,1].mean()),
        apex_parameter=float(sample_t[np.argmax(samples[:,1])]),apex_xy=samples[np.argmax(samples[:,1])].tolist()))

def run(input_path,out,smooth=1e-4,n_control=7,sensitivity=True):
    if not 5<=n_control<=9 or not 1e-6<=smooth<=1e-2:raise ValueError('Use 5-9 controls and smoothing between 1e-6 and 1e-2.')
    out=Path(out);out.mkdir(parents=True,exist_ok=True);raw=Path(input_path).read_bytes();data=json.loads(raw)
    obs,end_native,end_q,center,axes=observations(data);points=np.array([x['aligned'][:2] for x in obs]);ends=end_q[:,:2]
    fit=fit_arch(points,ends,n_control=n_control,smooth=smooth);curve=fit['curve']
    zref=float(end_q[:,2].mean());sample_t=np.linspace(0,1,501);xy=curve(sample_t);q=np.column_stack([xy,np.full(len(xy),zref)]);native=q@axes.T+center
    for i,o in enumerate(obs):
        foot=curve(fit['t'][i]);o.update(fitted_parameter=float(fit['t'][i]),projected_native=(np.r_[o['aligned'][:2],zref]@axes.T+center).tolist(),
            footpoint_xy=foot.tolist(),footpoint_native=(np.r_[foot,zref]@axes.T+center).tolist(),residual_native_units=float(fit['residuals'][i]),
            robust_weight=float(fit['robust_weights'][i]),effective_weight=float(o['weight']*fit['robust_weights'][i]))
    result=dict(algorithm_version='1.0',source_STL_sha256=data['source_sha256'],input_v2_json=str(Path(input_path).resolve()),input_v2_json_sha256=hashlib.sha256(raw).hexdigest(),
        method='Asymmetric endpoint-constrained robust cubic B-spline, with second-derivative regularization and Euclidean closest-point fitting',
        protocol='Current overall shape; 20 MD observations from P2-to-P2; projected M1 mesial points are hard endpoints. No symmetry constraint. No cusp/fossa/incisal-midpoint mixing.',
        clinical_status='Geometric reference candidate for subsequent space analysis; not an ideal treatment arch or a clinically validated available-space measurement.',
        units='Native STL coordinate units; physical millimetre calibration remains unconfirmed',
        input_tooth_ids=TOOTH_ORDER,excluded_from_soft_fit=['VL_M1','VR_M1','VL_M2','VR_M2'],
        pca_center=center.tolist(),pca_axes_columns=axes.tolist(),reference_plane_aligned_z=zref,
        reference_plane_definition='Existing V2 qXY plane, translated to mean qZ of the two M1 mesial boundary points. Endpoints and observations are orthogonally projected onto this plane.',
        endpoint_ids=['VL_M1_M','VR_M1_M'],endpoint_source_native=end_native.tolist(),endpoint_source_aligned=end_q.tolist(),
        endpoint_projected_native=(np.column_stack([ends,np.full(2,zref)])@axes.T+center).tolist(),
        control_points_xy=fit['control'].tolist(),control_points_native=(np.column_stack([fit['control'],np.full(n_control,zref)])@axes.T+center).tolist(),
        knots=fit['knots'].tolist(),degree=3,smoothing_lambda=smooth,n_control=n_control,
        huber_delta_native=fit['huber_delta_native'],all_teeth_equal_base_weight=True,
        projected_arch_length=fit['length'],length_definition='Integral of ||dC_xy/dt|| dt on [0,1]. The display curve is planar and has the same analytic 3-D length; source landmark height differences do not add length.',
        numerical_integration_error_estimate=fit['integration_error'],shape_diagnostics=fit['shape_diagnostics'],
        convergence=dict(converged=fit['converged'],iterations=fit['iterations'],optimizer=fit['optimizer'],warmup_history=fit['history']),
        observations=obs,curve_samples_native=native.tolist(),curve_samples_xy=xy.tolist(),sample_parameters=sample_t.tolist(),
        review_required=True,limitations=['Uses accepted candidate MD landmarks as proxies for the contact-measurement trajectory; these are not necessarily current physical neighbor contacts.',
            'Smoothing and curve family define the reference; residuals and parameter sensitivity are not anatomical error estimates.',
            'The existing V2 reference plane is retained; the clinical occlusal plane and ideal tooth arrangement are not inferred.',
            'Boundary points stay fixed; geometric checks must be reviewed independently for each fitted arch.',
            'No crowding or spacing diagnosis is made by this fitting stage.'])
    tests=[]
    if sensitivity:
        dense_t=np.linspace(0,1,10001);base_dense=curve(dense_t);base_tree=cKDTree(base_dense)
        omit_tooth=obs[int(np.argmax(fit['residuals']))]['tooth']
        for name,n,lam,exclude in [('lambda_half',n_control,smooth/2,None),('lambda_double',n_control,smooth*2,None),('six_controls',6,smooth,None),('eight_controls',8,smooth,None),('omit_largest_residual_tooth',n_control,smooth,omit_tooth)]:
            keep=np.array([o['tooth']!=exclude for o in obs]);f2=fit_arch(points[keep],ends,n_control=n,smooth=lam)
            variation=float(np.max(np.linalg.norm(f2['curve'](sample_t)-xy,axis=1)))
            other_dense=f2['curve'](dense_t)
            sampled_hausdorff=max(base_tree.query(other_dense)[0].max(),cKDTree(other_dense).query(base_dense)[0].max())
            tests.append(dict(name=name,n_control=n,smoothing_lambda=lam,omitted_tooth=exclude,length=f2['length'],length_change=f2['length']-fit['length'],
                max_same_parameter_curve_change=variation,max_bidirectional_sampled_curve_distance=float(sampled_hausdorff),distance_sampling_points=10001,
                passes_sampled_no_reverse_bending=f2['shape_diagnostics']['max_signed_curvature']<0,
                converged=f2['converged'],shape_diagnostics=f2['shape_diagnostics']))
    result['sensitivity']=tests
    (out/'reference_arch_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    with (out/'reference_arch_points.csv').open('w',encoding='utf-8-sig',newline='') as file:
        writer=csv.writer(file);writer.writerow(['point_index','t','x_native','y_native','z_native'])
        for i,(t,p) in enumerate(zip(sample_t,native)):writer.writerow([i,t,*p])
    with (out/'reference_arch_observations.csv').open('w',encoding='utf-8-sig',newline='') as file:
        writer=csv.writer(file);writer.writerow(['landmark','tooth','qX','qY','qZ_source','closest_parameter','normal_distance_native_units','robust_weight'])
        for o in obs:writer.writerow([o['id'],o['tooth'],*o['aligned'],o['fitted_parameter'],o['residual_native_units'],o['robust_weight']])
    with (out/'reference_arch_spline_controls.csv').open('w',encoding='utf-8-sig',newline='') as file:
        writer=csv.writer(file);writer.writerow(['control_index','x_native','y_native','z_native','qX','qY'])
        for i,(p,xy) in enumerate(zip(result['control_points_native'],result['control_points_xy'])):writer.writerow([i,*p,*xy])
    print(json.dumps(dict(length=fit['length'],converged=fit['converged'],iterations=fit['iterations'],RMS_residual=float(np.sqrt(np.mean(fit['residuals']**2))),max_residual=float(fit['residuals'].max()),shape=fit['shape_diagnostics'],sensitivity=tests),indent=2))
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input-v2',type=Path,required=True);parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--smooth',type=float,default=1e-4);parser.add_argument('--controls',type=int,default=7);parser.add_argument('--skip-sensitivity',action='store_true')
    args=parser.parse_args()
    if not 5<=args.controls<=9 or not 1e-6<=args.smooth<=1e-2:raise ValueError('Use 5–9 controls and smoothing between 1e-6 and 1e-2.')
    run(args.input_v2,args.output_dir,args.smooth,args.controls,not args.skip_sensitivity)
