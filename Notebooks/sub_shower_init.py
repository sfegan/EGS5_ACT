# coding: utf-8
base_rng = EGS5Simulations.EGS5RanluxSimpleRNG.instance()
pd_rng = EGS5Simulations.PredefinedDeviateSimpleRNG(base_rng);
pd_rng.setPD(0,mfpfrac);

detA = EGS5Simulations.PruningDetector(eprune)
detB = EGS5Simulations.EGS5SimpleIACTArray(atm, nlayer, emax, bfield, z0, ztop, nmedia)
detC = EGS5Simulations.TrackCountingDetector()
det = EGS5Simulations.EGS5UIDelegator()
det.appendDelegatee(detA)
det.appendDelegatee(detB)
det.appendDelegatee(detC)

layers = EGS5Simulations.VecLayer()
detB.getLayers(layers);
eff = EGS5Simulations.TelescopeEfficiency()
eff.scaleEffFromFile('Parameters/corsika_mirreff.dat')
eff.scaleEffFromFile('Parameters/corsika_quanteff.dat')
atmabs = EGS5Simulations.AtmosphericAbsorption('Parameters/corsika_atmabs.dat')
actyield = atmabs.integrateYield(z0, w0, eff)
for iscope in scopes:
    s = EGS5Simulations.EGS5ACTArrayImagingScope()
    s.x.set(iscope[0]*100, iscope[1]*100, z0)
    s.r   = iscope[2]*100 * 0.5
    s.res = res
    s.fov = fov
    s.setZnAz(zn/180.0*pi,0)
    s._yield = actyield
    detB.addScope(s);

egs5 = EGS5Simulations.EGS5System.instance(det,pd_rng)
egs5.initializeEGS5()
