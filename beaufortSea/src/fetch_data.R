read_ostracode <- function(csv_in){
  ostracode_raw <- read_excel(path = csv_in, 
                                  sheet = "T2 Ost MC29bin5G30bin10J32bin20",
                                  skip = 1, 
                                  range = cell_cols(c(7, 70:131)))
  # clean up column names
  ostracode_data <- ostracode_raw[,c(1,65:123)]
  names(ostracode_data) <- substr(names(ostracode_data), 1, 20)
  
  return(ostracode_data)
}

read_foram <- function(csv_in){
  foram_raw <- read_excel(path = csv_in, , 
                           sheet = "HLY1302 FORAM counts as valu",
                           range = cell_cols(c(7, 70:131)))
  
  # clean up column names
  foram_data <- foram_raw[c(1:25, 27:120),c(2,24:42)] 
  names(foram_data) <- substr(names(foram_data), 1, 20)
  
  return(foram_data)
}