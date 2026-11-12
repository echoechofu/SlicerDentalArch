"""Consistent P2-to-P2 required space versus M1-mesial reference arc."""
from pathlib import Path
import argparse,json,csv,hashlib
ORDER=['VL_P2','VL_P1','VL_C','VL_I2','VL_I1','VR_I1','VR_I2','VR_C','VR_P1','VR_P2']
def run(md_path,arch_path,out):
 md=json.loads(Path(md_path).read_text(encoding='utf8'));arch=json.loads(Path(arch_path).read_text(encoding='utf8'))
 if md['source_sha256']!=arch['source_STL_sha256']:raise ValueError('MD and arch data have different source STL digests.')
 if hashlib.sha256(Path(md_path).read_bytes()).hexdigest()!=arch['input_v2_json_sha256']:raise ValueError('Arch was not fitted from this exact MD result file.')
 byid={t['id']:t for t in md['teeth']};R=sum(byid[k]['width'] for k in ORDER);A=arch['projected_arch_length'];C=R-A
 pairs=[dict(tooth=k,width=byid[k]['width'],flags=byid[k]['review_flags']) for k in ORDER]
 perturb=[]
 for name in ['angle_-5_deg','angle_+5_deg','half_step','lower_height_+0.50']:
  if not all(name in byid[k]['sensitivity'] for k in ORDER):continue
  r=sum(byid[k]['sensitivity'][name]['width'] for k in ORDER)
  perturb.append(dict(name=name,required_space_R=r,change_R=r-R,C_with_fixed_baseline_A=r-A,scope='All10 widths changed coherently; A held fixed to isolate width-rule sensitivity, not an end-to-end clinical interval'))
 result=dict(source_sha256=md['source_sha256'],units='native STL coordinate units; no conversion to millimetres',
  protocol='Current-form asymmetric smooth reference arch between projected M1 mesial boundaries; sum intrinsic MD-v2 widths of the10 teeth P2-to-P2',
  range_tooth_ids=ORDER,boundary_landmarks=arch['endpoint_ids'],excluded_widths=['VL_M1','VR_M1','VL_M2','VR_M2'],
  required_space_R=R,available_space_A=A,signed_crowding_R_minus_A=C,space_discrepancy_A_minus_R=-C,
  positive_shortage=max(C,0.),positive_surplus=max(-C,0.),
  interpretation='Net space surplus under this geometric protocol' if C<0 else 'Net space shortage under this geometric protocol',
  formula='R=sum_i ||M_i-D_i||; A=integral_0^1 ||gamma_prime(t)||dt; signed_crowding C=R-A; space_discrepancy D=A-R',
  widths=pairs,width_rule_sensitivity_with_fixed_A=perturb,
  arch_parameter_sensitivity_with_fixed_R=[dict(name=t['name'],A=t['length'],C=R-t['length']) for t in arch['sensitivity']],
  limitations=['This is a signed net discrepancy for the chosen protocol, not the sum of individual tooth displacements or actual interdental gaps.',
   'A net surplus can coexist with local rotation or displacement. A small result is sensitive to landmark/plane/reference-curve choices.',
   'MD candidates approximate intrinsic crown width; not necessarily the current physical contacts between adjacent teeth.',
   'No ground-truth anatomical annotations or physical-scale calibration; no clinical severity grading.',
   'C is signed: negative means net surplus. Positive shortage is separately clamped to zero.'])
 out=Path(out);out.mkdir(parents=True,exist_ok=True)
 (out/'crowding_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
 with (out/'crowding_measurements.csv').open('w',encoding='utf-8-sig',newline='') as file:
  writer=csv.writer(file);writer.writerow(['quantity_or_tooth','native_coordinate_units','meaning'])
  for t in pairs:writer.writerow([t['tooth'],t['width'],'included intrinsic MD width'])
  for label,val,meaning in [('R',R,'required space'),('A',A,'reference arch available space'),('C=R-A',C,'positive shortage; negative surplus'),('D=A-R',-C,'positive surplus; negative shortage')]:writer.writerow([label,val,meaning])
 print(json.dumps(dict(R=R,A=A,C=C,D=-C,perturb=perturb),indent=2));return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--md',required=True,type=Path);p.add_argument('--arch',required=True,type=Path);p.add_argument('--output-dir',required=True,type=Path);a=p.parse_args();run(a.md,a.arch,a.output_dir)
