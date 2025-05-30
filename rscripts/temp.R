
makeOfName <- function(specDat)
{
  ret = specDat$metaData$parentMs2
  ret = paste(ret, "_", specDat$metaData$sequence,
              "_", specDat$metaData$scanNum, sep = "")
  if(specDat$metaData$precursorCharge != 0)
  {
    ret = paste(ret, "_", specDat$metaData$precursorCharge, sep = "")
  }
  ret = paste(ret, OF_EXT, sep = "")
  return(ret)
}
