mmhc_calc <- function(csv_location) {
  
  library(bnlearn)
  data = read.csv(csv_location, fileEncoding="UTF-8", sep = ",", dec = ",", header = FALSE, stringsAsFactors=FALSE) 
  data_num = as.data.frame(apply(data, 2, as.numeric))  #convertimos a variables numericas 
  
  
  dag_mmhc = mmhc(data_num, maximize.args = list(score = 'bic-g'))
  return (as.data.frame(dag_mmhc[["arcs"]]))
  
}




tabu_calc <- function(csv_location) {
  
  library(bnlearn)
  data = read.csv(csv_location, fileEncoding="UTF-8", sep = ",", dec = ",", header = FALSE, stringsAsFactors=FALSE) 
  data_num = as.data.frame(apply(data, 2, as.numeric))  #convertimos a variables numericas
  
  
  dag_tabu = tabu(data_num, score = 'bic-g')
  return (as.data.frame(dag_tabu[["arcs"]]))
  
}
