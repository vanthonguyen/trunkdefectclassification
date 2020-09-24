
The structure of the folder is:
data/
    treename/
            treename-id1-defect_type.xyz
            treename-id2-defect_type.xyz
            ...
centerline/
          treename/
                  centerline.xyz
                  radius


Each subfolder of data folder contains the defects of a tree (one file per defect). The defect type can be deduct from
file name:

et: paper and push pin used to identify defect type (see https://doi.org/10.1016/j.compag.2020.105332)
bra: sequential branch
gm: epicormic branch
cigm: epicormic branch merged with branch scars
ci: branch scar
br, br2, br3: burl
ama: bud cluster
sp: sphaeroblast
pi: picot
ec: bark
bru: noise


centerline folder contains centerline and radius of each centerline segment of tree corresponding in the data folder.

