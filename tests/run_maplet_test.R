########################################################################################
########################### --- maplet testing framework --- ###########################
########################################################################################

# SET VARIABLE
# script file in project test directory
file <- "/Users/jak2043/Box/results/Kelsey/test_framework/proj_tests_jan/stats_pipe_ex_test.R"

########################################################################################
################################### --- RUN TEST --- ###################################
########################################################################################
# load test and comparison objects
maplet:::mti_logstatus("Loading reference object...")
test_dir <- dirname(file)
obj_pre <- gsub(".R", "", basename(file))
test_file <- file.path(test_dir, paste0(obj_pre, "_test_steps.R"))
obj_dir <- file.path(test_dir, "REFERENCE_OBJECTS")
obj_files <- dir(obj_dir)
last_step <- stringr::str_extract(obj_files, "[[:digit:]]+") %>% gsub(".rds", "", .) %>% as.numeric %>% max()
test_object_file <- file.path(obj_dir, paste0(obj_pre, last_step, ".rds"))
D_REFERENCE <- readRDS(test_object_file)

maplet:::mti_logstatus("Running pipeline script...")
source(test_file)

# globally compare two SE objects
global_results_df <- mtt_equal(D1 = D_REFERENCE, D2 = D)

# compare reference SE objects step-by-step
if(any(global_results_df=="Fail")){

  if(!file.exists(test_file)) stop("Test step file does not exist!")
  stepwise_results_df <- mtt_compare_stepwise(file = test_file, ref_path = obj_dir, obj_pre = obj_pre)
  
}


