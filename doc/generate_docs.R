# generate Rd docs
devtools::document()

# convert RD into pdf
gen_rd_pdf <- "R CMD Rd2dvi --pdf --title='Test of foo' -o /tmp/foo.pdf man/*.Rd"
system(gen_rd_pdf)
