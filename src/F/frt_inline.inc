C
C This file contains functions which are inlined
C (when this filemin included) for better performance
C
      QFuni(x) = exp(-x)
      QFloc(x) = exp(-x)
C
C  Inverse functions
C
      QIuni(x) = -log(1.0e-35+x)
      QIloc(x) = -log(1.0e-35+x)
