### building the spatialdecon plugin

The basic idea: this folder contains the tools to build a DSPDA plugin running the spatialdecon package. 

The script "build spatialdecon plugin.R" assembles a single R script containing a main() function and
 all the code from the spatialdecon package. 
 
The script "testspatialdecon plugin.R" runs the assembled plugin. 

To test custom arguments, just save a copy of the plugin script, "SpatialDecon_plugin.txt", and modify the
 arguments up top as desired. 