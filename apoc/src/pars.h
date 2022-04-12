c=============================================================
c Parameters uses by P-align
c=============================================================
c
c PDB file related parameters
c
      integer        maxa,maxr,maxr2,maxs,maxc,maxk,maxp,maxt
      parameter    (
     -               maxa=100,    !max atoms per residue
     -               maxr=3000,   !max number of residues in PDB file
     -               maxr2=10000, !for display seq aln, double maxr
     -               maxs=100000, !max number of atoms in PDB file
     -               maxc=200,    !max number of contacts per residue
     -               maxk=10,     !max number of chains	
     -               maxp=100,    !max number of pockets/fragments per structure
     -               maxt=100000   !max number of structures
     -             )

c=============================================================
c
c Parameters for dynamic programming
c
      real   inf, ninf
      parameter    (
     -      	     inf=1.0E9,    !positive infinity
     -	             ninf=-1.0E9   !negative infinity
     -             )

c=============================================================
