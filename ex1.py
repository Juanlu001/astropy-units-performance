from astropy import units as u

from poliastro.twobody import Orbit
from poliastro.bodies import Earth
from poliastro.twobody.propagation import cowell

from poliastro.bodies import Body
from poliastro.twobody.rv import RVState

r0 = [-2384.46, 5729.01, 3050.46] * u.km
v0 = [-7.36138, -2.98997, 1.64354] * u.km / u.s

initial = Orbit.from_vectors(Earth, r0, v0)

@profile
def accel(t0, state, k):
    v_vec = state[3:]
    norm_v = (v_vec * v_vec).sum() ** .5
    return 1e-5 * v_vec / norm_v

@profile
def accel_slow(t0, state, k):
    r_vec, v_vec = state[:3], state[3:]
    _k = k * (u.km ** 3 / u.s ** 2)
    body = Body(None, _k, "_Dummy")
    _r = r_vec * u.km
    _v = v_vec * u.km / u.s
    ss = RVState(body, _r, _v)
    norm_v = (v_vec * v_vec).sum() ** .5
    return 1e-5 * v_vec / norm_v


km3s2 = u.km ** 3 / u.s ** 2
kms = u.km / u.s

@profile
def accel_so_so(t0, state, k):
    r_vec, v_vec = state[:3], state[3:]
    _k = k * km3s2
    body = Body(None, _k, "_Dummy")
    _r = r_vec * u.km
    _v = v_vec * kms
    ss = RVState(body, _r, _v)
    norm_v = (v_vec * v_vec).sum() ** .5
    return 1e-5 * v_vec / norm_v



@profile
def main():
    #initial.propagate(3 * u.day, method=cowell, ad=accel)
    #initial.propagate(3 * u.day, method=cowell, ad=accel_slow)
    initial.propagate(3 * u.day, method=cowell, ad=accel_so_so)


if __name__ == '__main__':
    main()
