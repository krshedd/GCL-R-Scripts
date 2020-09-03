TissueLocations_2R.GCL <- function(unit, username, password, bad_location = FALSE, all_data = FALSE) {
  #########################################
  # This function pulls tissue data from OceanAK, for a given storage location, and updates the tissue location maps.
  # Created by: Chase Jalbert
  # Created on: 9/1/2020
  #
  # Inputs~~~~~~~~~~~~~~~~~~
  #  unit - the storage unit(s) you're interested in; can be single or multiple units (e.g., "99" or c("WA","7","A"))
  #  username - your LOKI username and password
  #  password - your LOKI password
  #  bad_locations - TRUE or FALSE, do you want a list of tissues with incorrect location information; default = FALSE
  #  all_data - TRUE or FALSE, do you want a list of ALL tissues [potentially large file]]; default = FALSE
  #
  #
  # Outputs~~~~~~~~~~~~~~~~~
  #  tissuemap - an object mapping out all the tissues for a storage unit
  #  all_data - an object consisting of all tissues within OceanAK, within a storage unit
  #  bad_location - an object consisting of all tissues within OceanAK, containing unexpected, missing, or otherwise incorrect location information (i.e., SHELF_RACK and/or SLOT)
  #
  #
  # Example~~~~~~~~~~~~~~~~~
  #  username = "awesomeuser"
  #  .password = "awesomepassword1"
  #  unit = LETTERS[1:10] # we can cheat here, since B7 is A - Z... 
  #  bad_location = TRUE # yes, I want to identify all wrong locations
  #  all_data = TRUE # yes, I want all the data
  #
  #  TissueLocations_2R.GCL(unit = unit, username = username, password = .password, bad_location = bad_location, all_data = all_data)
  #  
  #  write_csv(x = tissuemap, path = paste0("C:/Users/csjalbert/Desktop/tissuemap_", Sys.Date(), ".csv")) # export CSV of the tissue map
  #  write_csv(x = all_data, path = paste0("C:/Users/csjalbert/Desktop/tissues_all_locations_", Sys.Date(), ".csv")) # export CSV of ALL tissues
  #  write_csv(x = bad_location, path = paste0("C:/Users/csjalbert/Desktop/bad_tissue_locations_", Sys.Date(), ".csv")) # export CSV of incorrect tissues
  ##########################################
  
  # Setup  
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(RJDBC, tidyverse, lubridate) #Install packages, if not in library and then load them.
  
  # Database setup
  ## This copies the "odbc8.jar" file to the R folder on your computer if it doesn't exist there. This file contains the java odbc drivers needed for RJDBC
  if (!file.exists(path.expand("~/R"))) {
    dir <- path.expand("~/R")
    
    dir.create(dir)
    
    bool <-
      file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))
    
  } else {
    if (!file.exists(path.expand("~/R/ojdbc8.jar"))) {
      bool <-
        file.copy(from = "V:/Analysis/R files/OJDBC_Jar/ojdbc8.jar", to = path.expand("~/R/ojdbc8.jar"))
      
    }
    
  }
  
  start.time <- Sys.time()
  
  options(java.parameters = "-Xmx10g")
  
  if (file.exists("C:/Program Files/R/RequiredLibraries/ojdbc8.jar")) {
    drv <-
      RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = "C:/Program Files/R/RequiredLibraries/ojdbc8.jar", " ")#https://blogs.oracle.com/R/entry/r_to_oracle_database_connectivity    C:/app/awbarclay/product/11.1.0/db_1/jdbc/lib
    
  } else {
    drv <-
      RJDBC::JDBC("oracle.jdbc.OracleDriver",
                  classPath = path.expand("~/R/ojdbc8.jar"),
                  " ")
    
  }
  
  ## Build database URL
  url <-
    LOKI_URL.GCL() # This is a function that gets the correct URL to access the database on the oracle cloud
  
  ## Connect to database
  con <-
    RJDBC::dbConnect(
      drv = drv,
      url = url,
      user = username,
      password = .password
    ) #The database connection
  
  ## Create the query
  gnoqry <-
    paste(
      "SELECT * FROM AKFINADM.V_GEN_SAMPLED_FISH_TISSUE WHERE UNIT IN (",
      paste0("'", unit, "'", collapse = ","),
      ")",
      sep = ""
    )
  
  # Data import
  ## Open the connection and pull data from the database
  dataAll <-
    RJDBC::dbGetQuery(conn = con, statement = gnoqry)
  
  ## Subset the data
  dataSubset <- dataAll %>%
    dplyr::filter(
      !grepl("ERICTEST", SILLY_CODE, ignore.case = TRUE),
      # drop Eric's test sillys
      PK_TISSUE_TYPE != "DNA",
      # drop DNA tissue type
      is.na(EXHAUSTED_HOW) |
        EXHAUSTED_HOW != "Discarded"
    ) %>% # drop tissues marked as discarded
    dplyr::select(
      c(
        FK_COLLECTION_ID,
        SILLY_CODE,
        FK_FISH_ID,
        STORAGE_ID,
        PK_TISSUE_TYPE,
        UNIT,
        SHELF_RACK,
        SLOT
      )
    ) %>% # select columns
    tidyr::unite(
      col = "tissue_id",
      SILLY_CODE,
      PK_TISSUE_TYPE,
      STORAGE_ID,
      sep = "_",
      remove = FALSE
    ) # create unique tissue identifier
  
  
  # Identify incorrect or missing tissue locations
  ## Note - SHELF_RACK and SLOT terminology differs with each storage location, so I've broken them out below
  bad_location <- dataSubset %>%
    tidyr::unite(col = "shelf_id", SHELF_RACK, SLOT, remove = FALSE) %>% # combine these into single column for an id
    # Find anything that does NOT match acceptable storage location code, for each unit
    dplyr::mutate(
      wrong = dplyur::case_when(
        # Freezer 99 - acceptable locations are: ###_UppercaseLetter
        UNIT %in% unit[stringr::str_detect(unit, "99")] &
          !stringr::str_detect(string = shelf_id, pattern = "[0-9]{3}_[A-Z]") ~ "wrong",
        # Warehouse - acceptable locations are: [1-25]_[1-10]
        UNIT %in% unit[stringr::str_detect(unit, "^(W[A-Z])$")] &
          !str_detect(string = shelf_id, pattern = "^([1-9]|1[0-9]|2[0-5])_([1-9]|10)$") ~ "wrong",
        # Freezer 7 or 8 - acceptable locations are: XXX?? UppercaseLetter_[1-9]
        UNIT %in% unit[stringr::str_detect(unit, "^([7-8])$")] &
          !stringr::str_detect(string = shelf_id, pattern = "^([A-Z])_([1-9])$") ~ "wrong",
        # B7 - acceptable locations are: XXX?? [-10]_UppercaseLetter
        UNIT %in% unit[stringr::str_detect(unit, "^([A-Z])$")] &
          !stringr::str_detect(string = shelf_id, pattern = "^([1-9]|10)_([A-Z])$") ~ "wrong",
        # Freezer 1 - acceptable locations are: XXX?? UppercaseLetter_[1-9]LowercaseLetter
        UNIT %in% unit[stringr::str_detect(unit, "1")] &
          !stringr::str_detect(string = shelf_id, pattern = "^([A-Z]_[1-9][a-z])$") ~ "wrong",
        TRUE ~ "correct"
      )
    ) %>%
    dplyr::filter(wrong == "wrong") %>%
    dplyr::select(-wrong)
  
  ## If, you said TRUE in setup, then assign the bad location object to your environment. From here you can export or do whatever you want with it.
  if (bad_location == TRUE) {
    assign(
      x = "bad_location",
      value = bad_location,
      pos = 1,
      envir = .GlobalEnv
    )
    message(paste(
      "All tissues with incorrect locations stored in object `bad_location`"
    ))
  }
  
  # Create final map
  ## Make fish number ranges
  fish_range <- dataSubset %>%
    dplyr::anti_join(
      bad_location,
      by = c(
        "FK_COLLECTION_ID",
        "tissue_id",
        "SILLY_CODE",
        "FK_FISH_ID",
        "STORAGE_ID",
        "PK_TISSUE_TYPE",
        "UNIT",
        "SHELF_RACK",
        "SLOT"
      )
    ) %>% # dropping any unknown/incorrect location collections
    dplyr::group_by(tissue_id) %>% # for each tissue id:
    dplyr::summarise(min = min(FK_FISH_ID), # find the min fish number
                     max = max(FK_FISH_ID)) # find the max fish number
  
  ## Assign the fish range to each fish
  TissueData <- fish_range %>%
    dplyr::left_join(
      dataSubset %>%
        dplyr::select(tissue_id, UNIT, SHELF_RACK, SLOT) %>% # Just need these columns for the map
        tidyr::unite(col = "shelf_id", SHELF_RACK, SLOT, remove = FALSE) %>% # Make a unique id column (same as above)
        # below filters for GOOD tissues, opposite of the bad_tissue section above
        dplyr::mutate(
          wrong = dplyur::case_when(
            # Freezer 99 - acceptable locations are: ###_UppercaseLetter
            UNIT %in% unit[stringr::str_detect(unit, "99")] &
              !stringr::str_detect(string = shelf_id, pattern = "[0-9]{3}_[A-Z]") ~ "wrong",
            # Warehouse - acceptable locations are: [1-25]_[1-10]
            UNIT %in% unit[stringr::str_detect(unit, "^(W[A-Z])$")] &
              !str_detect(string = shelf_id, pattern = "^([1-9]|1[0-9]|2[0-5])_([1-9]|10)$") ~ "wrong",
            # Freezer 7 or 8 - acceptable locations are: XXX?? UppercaseLetter_[1-9]
            UNIT %in% unit[stringr::str_detect(unit, "^([7-8])$")] &
              !stringr::str_detect(string = shelf_id, pattern = "^([A-Z])_([1-9])$") ~ "wrong",
            # B7 - acceptable locations are: XXX?? [-10]_UppercaseLetter
            UNIT %in% unit[stringr::str_detect(unit, "^([A-Z])$")] &
              !stringr::str_detect(string = shelf_id, pattern = "^([1-9]|10)_([A-Z])$") ~ "wrong",
            # Freezer 1 - acceptable locations are: XXX?? UppercaseLetter_[1-9]LowercaseLetter
            UNIT %in% unit[stringr::str_detect(unit, "1")] &
              !stringr::str_detect(string = shelf_id, pattern = "^([A-Z]_[1-9][a-z])$") ~ "wrong",
            TRUE ~ "correct"
          )
        ) %>%
        dplyr::filter(wrong == "wrong") %>%
        dplyr::select(-wrong),
      by = "tissue_id"
    ) %>%
    dplyr::distinct()
  
  ## Convert into map format
  data_map <- TissueData %>%
    tidyr::unite(range, min, max, sep = "-") %>%
    tidyr::separate(tissue_id,
                    into = c("silly", "tissue", "barcode"),
                    sep = "_") %>% # we have to split these to insert fish range
    tidyr::unite(map_id, silly, tissue, range, barcode, sep = " ") %>% # now making the final map entries
    dplyr::distinct() %>% # we just want unique entries
    dplyr::select(-c(shelf_id)) %>% # drop extra column
    dplyr::arrange(UNIT, SHELF_RACK, SLOT) %>% # sorting by shelf_rack, then slot
    dplyr::group_by(UNIT, SHELF_RACK, SLOT) %>%
    dplyr::group_modify( ~ {
      .x %>%
        mutate(row = row_number())
    }) # experimental magic... - modifies groups, based on previous group_by(); I am simply getting row numbers here
  
  ## Create tissue location map
  tissuemap <- data_map %>%
    tidyr::pivot_wider(names_from = SLOT,
                       values_from = map_id) %>%
    dplyr::select(-row) %>%
    dplyr::select(UNIT, SHELF_RACK, sort(tidyselect::peek_vars()))
  
  assign(
    x = "tissuemap",
    value = tissuemap,
    pos = 1,
    envir = .GlobalEnv
  )
  
  
  stop.time <- Sys.time()
  
  fulltime <- stop.time - start.time
  
  
  print(fulltime)
  message(paste0("Map of tissue locations stored in object 'tissuemap'"))
  
  
  if (all_data == TRUE) {
    #all_data %>% write_csv(path = paste0("tissues_all_locations_", Sys.Date(), ".csv")) # xxx do we want to just export a CSV to some directory?? desktop??
    assign(
      x = "all_data",
      value = dataAll,
      pos = 1,
      envir = .GlobalEnv
    )
    message(paste("All fish from database stored in object `all_data`"))
  }
}