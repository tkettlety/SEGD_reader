'''
`utils.py`: Utility functions for processing SEG-D data.
Includes:
- Function to rotate 3-component data (obspy stream) acquired using Sercel Galperin brackets.


MIT License

Copyright (c) 2024 Sacha Lapins

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''


import numpy as np
import obspy


def rotate_galperin_to_ENZ(stream, first_component_orientation="W", tilt_angle=35.26):
    """
    Rotate seismic data (obspy stream) from Galperin configuration (UVW) to ENZ.

    Convention is for first component (U) to face west and 2nd component (V) to be oriented 30 degrees clockwise from North.
    
    However, we deployed with one component oriented North, so I've included transformation for U oriented north and V oriented 120 degrees clockwise from north.

    Default tilt on each component is 35.26 degrees, as per original Galperin configuration.
    
    Returns:
    obspy stream with rotated data
    """

    alpha = np.pi*(tilt_angle / 180) # tilt angle from horizontal (default is 35.26 degrees)
    if first_component_orientation == 'N':
        # TAPIR configuration
        beta = np.pi*(120/180) # orientation of V (2nd) component wrt North (120 degrees clockwise)
        
        # Transformation matrix to convert to ENZ
        # I worked this out from scratch so might be wrong...
        transformation_matrix = np.array([
            [0, np.cos(alpha)*np.sin(beta), -np.cos(alpha)*np.sin(beta)],            # E
            [np.cos(alpha), np.cos(alpha)*np.cos(beta), np.cos(alpha)*np.cos(beta)], # N
            [np.sin(alpha), np.sin(alpha), np.sin(alpha)]                            # Z
        ])
    elif first_component_orientation == 'W':
        # Standard configuration
        beta = np.pi*(30/180) # orientation of V (2nd) component wrt North (30 degrees clockwise)

        # Transformation matrix to convert to ENZ
        # from https://en.wikipedia.org/wiki/Galperin_configuration
        transformation_matrix = np.array([
            [-np.cos(alpha), np.cos(alpha)*np.sin(beta), np.cos(alpha)*np.sin(beta)],  # E
            [0, np.cos(alpha)*np.cos(beta), -np.cos(alpha)*np.cos(beta)],              # N
            [np.sin(alpha), np.sin(alpha), np.sin(alpha)]                              # Z
        ])
    else:
        print("Unknown first component orientation. Only north ('N') and west ('W') currently supported. Doing nothing.")
        return None

    stations = list(set([tr.stats.station for tr in stream]))
    stations.sort()
    # Loop through stations and change those with three components for total duration of stream:
    for sta in stations:
        st = stream.select(station=sta)

        # Do checks (must have 3 components and non-zero data on all components)
        if (len(st) != 3) or (np.any([np.all(tr.data == 0) for tr in st])):
            continue

        st.sort() # Need component order to be sorted
        
        components_list = [tr.stats.component for tr in st]
        npts_list = [tr.stats.npts for tr in st]
        if (all(c in ['1', '2', '3', 'U', 'V', 'W'] for c in components_list)) & (all(n == npts_list[0] for n in npts_list)):
            # Try rotate data
            try:
                st.detrend("demean") # Mean removal probably necessary, as large offsets between sensors will overly scale combination?
                zne_data = np.dot(transformation_matrix, np.array([st[0].data, st[1].data, st[2].data]))
            except:
                continue
            
            # Remove traces from stream if successful
            for tr in st:
                stream.remove(tr)

            # Replace st data, change component name, then add back to st:
            st[0].data = zne_data[0,:]
            st[1].data = zne_data[1,:]
            st[2].data = zne_data[2,:]
            st[0].stats.channel = st[0].stats.channel[:-1] + 'E'
            st[1].stats.channel = st[1].stats.channel[:-1] + 'N'
            st[2].stats.channel = st[2].stats.channel[:-1] + 'Z'

            # Add traces back to stream
            for tr in st:
                stream.append(tr)

    stream.sort()
    return stream

