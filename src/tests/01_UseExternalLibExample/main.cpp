#include <iostream>
#include"MieosProjectOrganiserExtLib.hpp"

int main(int argc, char *argv[]){

   std::cout << std::endl << "MIEOSPO : Test 01 : Use of an external Library :" << std::endl << std::endl;

   //Call the static example function
   std::cout << "Test static function : " << std::endl;
   MieosProjectOrganiserExtLib::staticSayHi();
   std::cout << std::endl;

   //Call sayHi on the object
   std::cout << "Test sayHi function : " << std::endl;
   MieosProjectOrganiserExtLib* testObj = new MieosProjectOrganiserExtLib();
   testObj->sayHi();
   std::cout << std::endl;

   //Delete the object
   delete testObj;

   std::cout << "MIEOSPO : Test 01 : Ending.. " << std::endl << std::endl;

   return 0;
}
