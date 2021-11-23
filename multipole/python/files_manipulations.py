import  numpy as np


class pts:

    def __init__(self, p1, p2, p3, name):
        """Function to initialise a pts object.

        Args:
            p1 (numpy array): Array of the first coordinate for all the points.
            p2 (numpy array): Array of the second coordinate for all the points.
            p3 (numpy array): Array of the third coordinate for all the points.
            name (str): Name of the pts (without extension) for debugging
            purposes.
        """
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.size = p1.shape[0]
        self.name = name

    @classmethod
    def read(cls, pathname):
        """Function to generate a pts object from a file.

        Args:
            pathname (str): Full path including filename and extension.

        Returns:
            pts: pts object extracted from the file.
        """

        ptsfile = np.genfromtxt(pathname, delimiter=' ',
                                dtype=float, skip_header=1
                                )

        p1 = ptsfile[:, 0]
        p2 = ptsfile[:, 1]
        p3 = ptsfile[:, 2]

        full_name = pathname.split("/")[-1]
        name_noext_vec = full_name.split(".")
        name_noext = '.'.join(name_noext_vec[:-1])

        return cls(p1, p2, p3, name_noext)


class surf:

    def __init__(self, i1, i2, i3, tags=None):
        """Init function for the surf class.

        Args:
            i1 (numpy array of integers): Array with the first index for each
            element.
            i2 (numpy array of integers): Array with the second index for each
            element.
            i3 (numpy array of integers): Array with the third index for each
            element.
            mesh_from (str): Mesh where the vtx comes from, for debugging
            purposes.
            tags (numpy array of integers, optional): In case of being a
            surfmesh, the list of tag for each element. If it is a surf, takes
            the value None. Defaults to None.
        """
        self.i1 = i1.astype(int)
        self.i2 = i2.astype(int)
        self.i3 = i3.astype(int)
        if (tags is None):
            self.tags = tags
        else:
            self.tags = tags.astype(int)

        self.size = i1.shape[0]

    @classmethod
    def read(cls, pathname):
        """Function to read a .surf, .elem or .surfmesh file and convert it to a
        surf object.

        Args:
            pathname (str): Full path (including filename and extension).
            mesh_from (str): Mesh where the vtx comes from, for debugging
            purposes.
        Returns:
            surf: surf object extracted from the file.
        """

        surfmesh_str = pathname.split(".")[-1]
        if surfmesh_str == "surfmesh" or pathname.split(".")[-2][0:4] == "part":
            num_cols = (1, 2, 3, 4)
            is_surfmesh = True
        else:
            num_cols = (1, 2, 3)
            is_surfmesh = False

        surffile = np.genfromtxt(pathname, delimiter=' ', dtype=int, skip_header=True, usecols=num_cols)

        i1 = surffile[:, 0]
        i2 = surffile[:, 1]
        i3 = surffile[:, 2]

        if is_surfmesh:
            tags = surffile[:, 3]
        else:
            tags = None

        return cls(i1, i2, i3, tags)
#
#     def tovtx(self):
#         """Function to extract the vtx from a surf object, removing duplicates.
#
#         Returns:
#             vtx: vtx object with the indices from the surf.
#         """
#         vtxvec = np.unique(np.concatenate((self.i1, self.i2, self.i3),
#                                           axis=0))
#
#         return vtx(vtxvec, self.mesh_from)
#
#     def write(self,pathname):
#         """Function to write a surf object to a file.
#
#         Args:
#             pathname (str): Full path (including filename and extension).
#         """
#
#         header = np.array([str(self.size)])
#         elemtype = np.repeat("Tr",self.size)
#         data = [elemtype, self.i1, self.i2, self.i3]
#
#         if(self.tags is not None):
#             data = [elemtype, self.i1, self.i2, self.i3, self.tags]
#
#         np.savetxt(pathname, header, fmt='%s')
#
#         with open(pathname, "ab") as f:
#             np.savetxt(f, np.transpose(data), fmt = "%s")
#
#
# class vtx:
#
#     def __init__(self, indices, mesh_from):
#         """Init function for the vtx class.
#
#         Args:
#             indices (numpy array of integers): Array of indices.
#             mesh_from (str): Mesh where the vtx comes from, for debugging
#             purposes.
#         """
#         self.indices = indices
#         self.size = len(indices)
#         self.mesh_from = mesh_from
#
#     def write(self,pathname):
#         """Function to write a vtx object to a file.
#
#         Args:
#             pathname (str): Full path (including filename and extension).
#         """
#         array = np.array(self.indices)
#         header = [self.size, "intra"]
#         str_to_write = np.append(header, array)
#
#         with open(pathname, 'w') as f:
#             for item in str_to_write:
#                 f.write("%s\n" % item)