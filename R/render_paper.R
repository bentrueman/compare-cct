
library("tidyverse")

# paper:

path <- "Rmarkdown/paper.Rmd"
path_temp <- str_replace(path,"(.+/)(.+)", "\\1_\\2")

read_file(path) %>% 
  # move arabic numberal refs after punctuation:
  str_replace_all("\\s+(\\[@.*?\\])([.,;:])?", "\\2\\1") %>% 
  write_file(path_temp)

rmarkdown::render(path_temp)

file.rename(str_replace(path_temp, "\\.Rmd", "\\.docx"), str_replace(path, "\\.Rmd", "\\.docx"))

unlink(path_temp)

# supplement:

rmarkdown::render("Rmarkdown/paper-si.Rmd")
