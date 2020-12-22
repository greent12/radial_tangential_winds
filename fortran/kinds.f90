MODULE kinds

!Define type kinds in this module

IMPLICIT NONE
SAVE
INTEGER,PARAMETER :: single=selected_real_kind(p=6,r=37)
INTEGER,PARAMETER :: double=selected_real_kind(p=13)

END MODULE  
