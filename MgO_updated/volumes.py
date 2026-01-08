import numpy as np
volumes_m3gnet1 = np.arange(10.5,16.0,0.5)
volumes_m3gnet2 = np.arange(16.0, 21.7, 0.3)
volumes_m3gnet = np.concatenate((volumes_m3gnet1, volumes_m3gnet2))
volumes_m3gnet = np.append(volumes_m3gnet,19.190)
volumes_m3gnet = np.sort(volumes_m3gnet)
