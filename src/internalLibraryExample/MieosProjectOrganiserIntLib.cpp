#include "MieosProjectOrganiserIntLib.hpp"
#include <iostream>

MieosProjectOrganiserIntLib::MieosProjectOrganiserIntLib(){

}

MieosProjectOrganiserIntLib::~MieosProjectOrganiserIntLib(){

}

void MieosProjectOrganiserIntLib::staticSayHi(){
   std::cout <<"Static Hi from internal example class"<<std::endl;
}

void MieosProjectOrganiserIntLib::sayHi(){ 
   std::cout <<"Hi from internal example class"<<std::endl;
}
