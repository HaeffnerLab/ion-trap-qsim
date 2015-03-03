        #!/usr/bin/env python
from __future__ import division, absolute_import, print_function, unicode_literals
class Laser:
        
        def __init__(self, ion_num, freq, intensity):
                self.ion_num   = ion_num
                self.freq         = freq
                self.intensity    = intensity
                
                
class Pulse:

        def __init__(self, **kwargs):
                self.ion_num = kwargs['ion_num']
                self.sideband_num = kwargs['sideband_num']
                self.intensity = kwargs['intensity']
                self.duration  = kwargs['duration']

                self.eta = 0.05

