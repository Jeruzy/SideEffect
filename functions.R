#Function to call pubchem API and obtain inchi in txt format

get_inchi <- function(pid){
  
  url_base <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
  url_output <- "/property/InChI/txt"
  final_url <- paste(url_base, pid, url_output, sep = "")
  
  result <- httr::GET(final_url)
  content <- httr::content(result, "text")
  return(content)
  
}

vec_ptype2.glue <- function(x, y, ...) UseMethod("vec_ptype2.glue", y)
vec_ptype2.glue.character <- function(x, y, ...) x
vec_ptype2.glue.default <- function(x, y, ...) vctrs::vec_default_ptype2(x, y)
vec_ptype2.character.glue <- function(x, y, ...) y

vec_cast.glue <- function(x, to, ...) UseMethod("vec_cast.glue")
vec_cast.glue.character <- function(x, to, ...) glue::as_glue(x)
vec_cast.character.glue <- function(x, to, ...) unclass(x)