## Table of Contents
[Essential Packages For This Guide](#)
[Getting Data Into R](#)
[Functions For Data Exploration]()
[Data Transformation](#Data-Transformation-in-R)

## Essential Packages For This Guide:

>**Base R**
> **data.table** – fast, powerful, memory efficient, and succinct package for working with tabular data in R. Because it’s so succinct, it has a significant learning curve–see **[Intro to data.table Video](https://www.youtube.com/watch?v=uueVddWwbkk)** and [data.table Operators Video](https://www.youtube.com/watch?v=x6ufCO7H9pY) for a gentle introduction to data.table and use this **[data.table Cheat Sheet](https://res.cloudinary.com/dyd911kmh/image/upload/v1653830846/Marketing/Blog/data_table_cheat_sheet.pdf)** and the help documentation available through Rstudio for further reference while you code.
> **ggplot2**  

In this guide, functions that aren't a part of base R will usually be written with their package specified using " :: ", ex.  igraph::**distances**( ) . Multiple packages will sometimes use some functions of the same name, so this notation is used in R code to specify which package the desired function comes from--specifying the package is not usually necessary, but the notation will be used in this guide for clarity. Some packages with unusual syntax will not be written with this package-first notation (ex. data.table, which makes extensive use of square brackets).   


## Getting Data Into R:

* data.table::**fread**( ) and data.table::**fwrite**( ) 
	- These are usually the most time and memory efficient way to import and export tabular data from R. They are at least an order of magnitude faster than base r read and write functions. Not limited to data table objects, they can write matrices, data frames, and data tables.
	-  The file format of data.table::**fwrite**( ) is conveniently specified by the file extension in the filename you provide. Example: you want to save a data table object, shoeSizes, as a .tsv, you would write  
		 ```data.table::fwrite( shoeSizes, "shoe_sizes_table.tsv" ) ```  
	- These can read and write gzipped files without an additional compression or unpacking function to save storage space ( .gz extension ). Append ".gz" to the file extension you want to use, example: "table.csv.gz"

- vroom::**vroom**( ) 
	* This can handle extremely large files better, depending on how you need to use the information inside of them, since it only reads the content of the file as you access it. For example, if you have a .tsv with 1000 columns and a file size of 10GB, but you only need to access 4 columns, vroom would be ideal. 

* utils::**read.table**( ), utils::**read.csv**( ), or readxl::**read_excel**( )
	* These functions are slower and use more memory than **fread**( ), but they will return a data frame instead of a data table if you prefer that, and and the file size is not limiting. 
	* Additionally, readxl::**read_exce**l( ) can be very convenient for excel files, although due to the wide variety of excel file formats with the same extension, it does not always work. If you can access a plain delimited text file, that's ideal.
	* utils::**write.table**( ) and utils::**write.csv**( ) are their counterparts for writing files, but are almost always inferior to data.table::**fwrite**( ).

- **save**( ), **load**( ), and **readRDS**( )
	* These are for saving or importing R objects with .RData and .rds files, respectively. 
	* .RData files can contain multiple R objects, and .rds files contain only one object. For this reason, **load**( ) creates these objects in the current (or specified) R environment without an assignment operator while **readRDS**( ) returns an object that should be assigned to a variable.
	* BE CAREFUL, load will overwrite objects in your environment with the same name. 
	* Most useful for saving/importing complicated list objects that aren't easily formatted as tables, especially those used internally by packages like igraph. 

## R Functions for Easy Data Exploration:
Essential tools for quickly checking something in the console or for debugging. Text in *italics* describes the type of object and/or the content of the object expected. 

* utils::**str**( *dataframe or data.table* )
	- “Structure,” prints column names and first few items of each column on a new line. Gets messy fast with nested lists. Essential for quickly finding column names and column data formatting when you’re working with many tables. 

- utils::**head**( *dataframe, or data.table* )  or  utils::**tail**( ... )
	- Similar function to utils::**str**( ), but prints the first 6 rows or last 6 rows of a table to the console in a tabular structure. Gets messy fast with lots of columns. 

* **unique** ( *character or numeric vector* ) 
	* Show the unique categorical values in a vector
	- **length**( **unique**( *vector* ) ) can be used to count the number of unique categorical values in a vector. Can also be typed using pipes |>
		```\<vector> |> unique() |> length()```

-  *data.table*\[ , **.N**, by = .( *column name(s)*) ]
	- “.N” in the second argument/j argument of a data.table call will count the number of rows for each unique value of the “by” column(s) (or combination of values if multiple columns) and print the table.
	- It also returns a table you can assign, or you can use the “:=” assignment operator to add the column to the original data table, ex:
		 ``` data.table[ , count := .N , by = .( column name(s) ) ] ```

* utils::**View**( *function name* )
	- Open the source code of a function definition in R. Very useful for debugging or otherwise quickly navigating large bodies of code.

### Quick Plotting Functions

* graphics::**hist**( *numeric vector* )
	- “Histogram” plots a histogram of a numeric vector.
	- The argument breaks = *number of bins* is useful for increasing or decreasing the resolution of the histogram bins, and the argument xlim = c( *lower x lim*, *upper x lim*) is very useful for getting a better resolution on the range you’re interested in. More functionality can be found on the help page.   

- ggvenn::**ggvenn**( *named list of categorical vectors* )
	- gg "Venn diagram" plots the overlap between 2 or more categorical variables; the overlap is not to scale, the plot is labelled with the percentage and absolute number of shared and exclusive values. Names of the vectors in the list label the circles of the diagram. 



## Data Transformation in R
There are 2 main R packages for data transformation. nOe is the tidyverse package, which is a suite of particularly human-readable packages, including dplyr, tidyr, stringr, and forcats. The other is data.table. 

### Subsetting

### Assignment
#### Screening

### Merging

### Shaping

### Summarization



