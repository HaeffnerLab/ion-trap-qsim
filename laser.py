        #!/usr/bin/env python

from __future__ import division, absolute_import, print_function, unicode_literals
from numpy import pi

class Laser:
        
        def __init__(self, **kwargs):
                self.ion_num = kwargs['ion_num']
                self.sideband_num = kwargs['sideband_num']
                self.intensity =  kwargs['intensity']
                self.detuning  =  kwargs['detuning']
                self.phase     =  kwargs['phase']
                self.eta = 0.05


class Pulse:

        def __init__(self, **kwargs):
                self.ion_num = kwargs['ion_num']
                self.sideband_num = kwargs['sideband_num']
                self.intensity = kwargs['intensity']
                self.duration  = kwargs['duration']
                self.phase     =  kwargs['phase']
                self.detuning  =  kwargs['detuning']
                self.eta = 0.05

