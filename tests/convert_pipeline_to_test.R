########################################################################################
###################### --- setup for maplet testing framework --- ######################
########################################################################################

# SET VARIABLE
# maplet pipeline file to convert to test file; directory will be used as project test directory
file <- "/Users/jak2043/Box/results/Kelsey/test_framework/proj_tests_jan/stats_pipe_ex_test.R"

########################################################################################
################################### --- RUN SETUP --- ##################################
########################################################################################

library(magrittr)

# check if files exist
test_dir <- dirname(file)
obj_pre <- gsub(".R", "", basename(file))
test_step_file <- file.path(test_dir, paste0(obj_pre, "_test_steps.R"))
save_step_file <- file.path(test_dir, paste0(obj_pre, "_save_steps.R"))
if(!file.exists(file)) stop("Input file does not exist!")
if(file.exists(test_step_file)) stop("Test step file already exists!")
if(file.exists(save_step_file)) stop("Save step file already exists!")

# sub directory for reference objects
obj_dir <- file.path(test_dir, "REFERENCE_OBJECTS/")

code <- parse(file)
tokens <- as.list(code)

step_res_lst <- list()

step_res_idx <- 1
contains_mt <- FALSE
mt_count <- 1

for(ii in seq_along(tokens)){
  
  part <- tokens[[ii]]
  for(i in seq_along(part)){
    p <- as.character(part[i])
    # if no mt_ functions in line of code, rewrite to test script as is
    if(grepl("mt_", p)){
      # if multiple mt calls, split into individual calls
      if(grepl("%>%", p)){
        mt_calls <- strsplit(p, "%>%") %>% unlist()
        for(mtc in mt_calls){
          mtc <- stringr::str_trim(mtc)
          if(startsWith(mtc, "mt_")){
            if(mt_count > 1){
              step_res_lst[step_res_idx] <- paste0("D %<>% ", mtc)
              step_res_idx = step_res_idx + 1
            }else{
              step_res_lst[step_res_idx] <- paste0("D <- ", mtc)
              step_res_idx = step_res_idx + 1
            }
            # if a mt_plots function, remove results - plots are not compared and objects with plots are large
            if(startsWith(mtc, "mt_plots")){
              step_res_lst[step_res_idx] <- "metadata(D)$results[[length(metadata(D)$results)]]$output <- NULL"
              step_res_idx = step_res_idx + 1
            }
            step_res_lst[step_res_idx] <- glue::glue("saveRDS(D, file=\"{obj_dir}{obj_pre}{mt_count}.rds\")")
            step_res_idx = step_res_idx + 1
            mt_count = mt_count + 1
          }
        }
      }
      contains_mt <- TRUE
    } 
    
  }
  if(contains_mt==FALSE){
    p <- deparse(part)
    if(length(p) > 1) p <- paste0(p, collapse = "")
    step_res_lst[[step_res_idx]] <- p
    step_res_idx = step_res_idx + 1
  }
  
  
}

# get code without save commands
step_res_lst <- unlist(step_res_lst)
mt_funs_only <- sapply(step_res_lst, function(x){if(startsWith(x, "saveRDS")){NULL}else{x}}) %>% unlist() %>%
  sapply(., function(x){if(startsWith(x, "metadata")){NULL}else{x}}) %>% unlist()

# crash if last call is not mt_clean_remove_results(remove="plots")
clean_last = grepl("mt_clean_remove_results", mt_funs_only[length(mt_funs_only)]) &&  grepl("plots", mt_funs_only[length(mt_funs_only)])
if(!clean_last){
  stop("The command mt_clean_remove_results(remove=\"plots\") must be the last command in the pipeline.")
}

# create sub directory
if(dir.exists(obj_dir)){
  stop(glue::glue("Directory already exists: {obj_dir}"))
}else{
  dir.create(obj_dir)
}

# write save steps file
writeLines(step_res_lst, save_step_file, sep="\n")

# get code without save commands
# write test steps file
writeLines(mt_funs_only, test_step_file, sep="\n")

# source save step file
source(save_step_file)


