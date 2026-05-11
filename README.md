# CorsikaParallelism.jl                                                                                                                                             
                                                                                                                                                                      
  A Julia package for parallelizing [CORSIKA 8](https://www.iap.kit.edu/corsika/)                                                                                     
  air shower simulations. The primary shower is run until its first hadronic                                                                                          
  interaction, the resulting secondaries are extracted, and each is re-run as an                                                                                      
  independent sub-shower dispatched across multiple parallel workers.                                                                                                 
                                                                                                                                                                      
  ## How it works                                                                                                                                                     
                                                                                                                                                                      
  1. The primary particle is run through a modified CORSIKA 8 binary that stops
     after the first interaction and writes secondaries to a parquet file.                                                                                            
  2. The orchestrator reads the secondaries and distributes them as independent                                                                                       
     sub-shower jobs across `W` workers.                                                                                                                              
  3. Each sub-shower runs to observation level and its final-state particles are                                                                                      
     collected.                                                                                                                                                       
                                                                                                                                                                      
  This is a **1-generation recursion** scheme — sub-showers are not split further.
                                                                                                                                                                      
  ## Schedulers                                                                                                                                                       
                                                                                                                                                                      
  Four scheduling strategies are implemented in `src/schedulers.jl`:                                                                                                  
  - `naive` — assigns sub-showers in a round-robin manner                                                                                                                         
  - `threshold` — groups sub-showers by energy and assigns energy groups to workers                                                                                           
  - `binpack` — bin-packs sub-showers by estimated runtime to balance worker load                                                                                                
  - `workstealing` — workers pull from a shared queue; best overall performance                                                                                       
                                                                                                                                                                    
  ## Requirements                                                                                                                                                     
                                                                                                                                                                    
  - Julia 1.12+                                                                                                                                                       
  - A modified CORSIKA 8 binary with first-interaction stopping support                                                                                               
  - FLUKA (for hadronic interactions); set `FLUPRO` environment variable                                                                                              
  - `LD_LIBRARY_PATH` pointing to Julia's lib directory                                                                                                               
                                                                                                                                                                                                                                                                                                                    
  Results summary                                                                                                                                                     
                                                                                   
  Tested at 10 TeV, 100 TeV, and 1 PeV primary proton energies on the Harvard
  FASRC cluster. Workstealing with W=8 workers gives the best performance. Strong
  scaling holds up to ~8 workers and plateaus beyond that due to process-launch                                                                                       
  and I/O overhead. Final-state particle counts, total energy, and energy                                                                                             
  distributions agree with the serial CORSIKA 8 baseline within statistical                                                                                           
  uncertainties (KS statistic ~5–10%). 
