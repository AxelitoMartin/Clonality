get.oneseg.bdry <-
function(alpha, nperm) {
  max.ones <- floor(nperm*alpha) + 1
  sbdry <- DNAcopy::getbdry(eta=0.05, nperm, max.ones)
}

