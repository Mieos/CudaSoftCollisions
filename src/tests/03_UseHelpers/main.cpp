#include <iostream>
#include <string>
#include <fstream>

#include "pathsHelpers.hpp"

int main(int argc, char *argv[]){

   std::cout << std::endl << "MIEOSPO : Test 03 : Use of the helper :" << std::endl << std::endl;

   //Get and print the path to data
   std::cout << "Get the path to data. " << std::endl;
   std::string path = MIEOS_HELPERS_DATAPATH_GENERATED;
   path = path + "/dataExampleFolder/dataExampleFile.txt";
   std::cout << "The absolute path to data is : " << path << std::endl << std::endl;

   //Read a data file
   std::cout << "Reading the file : " << std::endl;
   std::ifstream mfile(path);
   std::string line;
   if (mfile.is_open()){
      
      while ( std::getline(mfile,line) ){
      
         std::cout << line << std::endl;
      }

   }
   std::cout << std::endl;

   std::cout << "MIEOSPO : Test 03 : Ending.. " << std::endl << std::endl;

   return 0;
}
