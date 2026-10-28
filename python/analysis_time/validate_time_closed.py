import json
from pathlib import Path
import numpy as np
from scipy.optimize import nnls
from numpy.linalg import norm,lstsq,svd

# Exact differentiation matrix of the periodic trigonometric interpolant,
# odd N, dimensionless tau=omega*t. No Fourier transform is calculated.
N=31
H=(N-1)//2
tau=2*np.pi*np.arange(N)/N
jj=np.arange(N)[:,None]-np.arange(N)[None,:]
D=np.zeros((N,N)); off=jj!=0
D[off]=.5*((-1.)**jj[off])/np.sin(np.pi*jj[off]/N)
# Sample-defined nonsinusoidal voltage. The saw component keeps all finite
# modes appreciably excited; the exact claim concerns this finite interpolant.
saw=np.mod(tau+0.37,2*np.pi)/np.pi-1
raw=np.tanh(1.6*np.sin(tau))/np.tanh(1.6)+0.055*saw
raw-=raw.mean()
u=230*raw/np.sqrt(np.mean(raw**2))
omega=2*np.pi*50
Rload=7.
Lload=.0167
i_rl=np.linalg.solve(Rload*np.eye(N)+Lload*omega*D,u)
u_shift=np.roll(u,2)/max(abs(u))
i_nonlinear=5*np.sign(u_shift)*abs(u_shift)**5
i=i_rl+i_nonlinear
i-=i.mean()

# Full sample model rates are 1,...,H. Pole selection uses the known
# differentiation operator, not an input Fourier spectrum.
lam=np.arange(1,H+1,dtype=float)**2
nu=np.arange(1,H,dtype=float)+.5
mu=nu**2
# Columns have common units volts; positive coefficients have units siemens.
phi_C=D@u
phi_L=lstsq(D,u,rcond=None)[0]
Phi=np.column_stack([phi_C,phi_L]+[np.linalg.solve(D@D+m*np.eye(N),D@u) for m in mu])
# Y0(z)=prod(z^2+lam)/(z*prod(z^2+mu)); its positive Foster residues
# span the unique positive null vector of the H by H+1 response map.
zC=1.
zL=float(np.prod(lam)/np.prod(mu))
zres=[]
for k,m in enumerate(mu):
    zres.append(float(np.prod(lam-m)/(-m*np.prod(np.delete(mu,k)-m))))
z=np.r_[zC,zL,zres]
assert np.all(z>0)
a0=lstsq(Phi,-i,rcond=None)[0]
alpha=float(np.max(-a0/z))
a=a0+alpha*z
# Round-off clipping only, not a KKT substitute.
a[np.abs(a)<1e-12*max(abs(a))]=0
q=Phi@a
# Independent physical-current nonnegative least squares.
a_nnls,res_nnls=nnls(Phi,-i,maxiter=10000)
q_nnls=Phi@a_nnls
# Orthogonal projection using a direct SVD basis, not derivative powers.
U,s,Vh=svd(Phi,full_matrices=False)
rank=int(np.count_nonzero(s>s[0]*1e-11))
Q=U[:,:rank]
q_orth=-Q@(Q.T@i)

rms=lambda x:float(np.sqrt(np.mean(x*x)))
power=lambda x:float(np.mean(u*x))
P=power(i)
ia=P/np.mean(u*u)*u
# SI component values for normalized coefficients.
components={"C_shunt_F":float(a[0]/omega),"L_shunt_H":float(1/(omega*a[1])) if a[1]>0 else None,"branches":[]}
for k,ak in enumerate(a[2:]):
    components["branches"].append({"pole_hz":float(50*nu[k]),"coefficient_S":float(ak),"L_H":float(1/(omega*ak)) if ak>0 else None,"C_F":float(ak/(omega*mu[k])) if ak>0 else None})
result={
"N":N,"rank":rank,"column_count":Phi.shape[1],"omega_rad_s":omega,
"differentiator_skew_relative":float(norm(D+D.T)/norm(D)),
"differentiator_sine_derivative_max_error":float(max(abs(D@np.sin(tau)-np.cos(tau)))),
"voltage_rms_V":rms(u),"power_W":P,"current_before_A":rms(i),
"current_after_closed_A":rms(i+q),"current_after_nnls_A":rms(i+q_nnls),
"current_after_projection_A":rms(i+q_orth),"fryze_bound_A":rms(ia),
"compensator_power_W":power(q),
"power_delta_W":power(i+q)-P,
"positive_null_min":float(min(z)),"positive_null_max":float(max(z)),
"positive_null_relative_residual":float(norm(Phi@z)/(norm(Phi)*norm(z))),
"alpha":alpha,"coefficient_min":float(min(a)),"coefficient_max":float(max(a)),
"active_components":int(np.count_nonzero(a>1e-12*max(a))),
"projection_current_max_error_A":float(max(abs(q-q_orth))),
"nnls_current_max_error_A":float(max(abs(q-q_nnls))),
"nonzero_singular_min":float(s[rank-1]),"singular_max":float(s[0]),
"component_values":components,
"notes":["Exact finite-dimensional periodic trigonometric interpolation model; no FFT used.","Voltage is defined directly by nonlinear/sample operations. All poles are chosen from the known derivative operator's interlacing grid.","The positive-null shift proves passive realizability in this ideal model, without component stress/tolerance constraints.","Component values can be impractical; this calculation does not establish robustness of an analog design."]
}
Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2))
np.savez(Path(__file__).with_name('validation_arrays.npz'),tau=tau,u=u,i=i,D=D,Phi=Phi,z=z,a0=a0,a=a,a_nnls=a_nnls,q=q,q_orth=q_orth,q_nnls=q_nnls)
print(json.dumps({k:v for k,v in result.items() if k not in ('component_values','notes')},indent=2))
