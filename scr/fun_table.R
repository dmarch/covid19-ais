#-----------------------------------------------------------------------
# fun_table         Suite of function to generate a docx table
#-----------------------------------------------------------------------


# create new word document
new.word.doc=function(){
  my.doc=read_docx()
  return(my.doc)
}

# add an empty line
add.empty.line=function(doc){
  body_add_par(doc, " ")
  return("empty line added")
}

#add page break
add.page.break=function(doc){
  body_add_break(doc, pos="after")
  return("page break added")
}

# start landscape
start.landscape=function(doc){
  doc=body_end_section_continuous(doc)
  return("landscape orientation started")
}


# end landscape
end.landscape=function(doc){
  doc=body_end_section_landscape(doc)
  return("landscape orientation ended")
}

# add a title
add.title=function(doc, my.title){
  my.prop=fp_text(font.size = 12, bold = TRUE, font.family = "Arial")
  the.title=fpar(ftext(my.title, prop=my.prop))
  body_add_fpar(doc, the.title)
  #body_add_par(doc, " ")
  return("title added")
}

# add a text
add.text=function(doc, my.text){
  my.prop=fp_text(font.size = 10, bold = FALSE, font.family = "Arial")
  the.title=fpar(ftext(my.text, prop=my.prop))
  body_add_fpar(doc, the.title)
  body_add_par(doc, " ")
  return("text added")
}

# add an image, such as a jpg or png file
add.image=function(doc, image, h=5, w=5){
  body_add_img(doc, src=image,
               height=h, width=w,
               style="centered")
  return("image added")
}


# add a data frame as a table
add.table=function(doc, tbl){
  # add the table to the document
  flextable::body_add_flextable(doc, 
                                value = tbl, 
                                align = "left" )
  return("table added")
}
