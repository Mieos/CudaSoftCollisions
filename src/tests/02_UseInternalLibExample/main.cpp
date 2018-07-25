#include <iostream>
#include"MieosProjectOrganiserIntLib.hpp"

int main(int argc, char *argv[]){

   std::cout << std::endl << "MIEOSPO : Test 02 : Use of an internal Library :" << std::endl << std::endl;

   //Call the static example function
   std::cout << "Test static function : " << std::endl;
   MieosProjectOrganiserIntLib::staticSayHi();
   std::cout << std::endl;

   //Call sayHi on the object
   std::cout << "Test sayHi function : " << std::endl;
   MieosProjectOrganiserIntLib* testObj = new MieosProjectOrganiserIntLib();
   testObj->sayHi();
   std::cout << std::endl;

   //Delete the object
   delete testObj;

   std::cout << "MIEOSPO : Test 02 : Ending.. " << std::endl << std::endl;

   return 0;
}
