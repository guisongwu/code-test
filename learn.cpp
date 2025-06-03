// ==================================== 命令行参数 =======================================
/* #include <iostream> */
/* using namespace std; */

/* int main(int argc, char* argv[]) { */
/*     cout << "参数数量: " << argc - 1 << endl; */
/*     cout << "程序名: " << argv[0] << endl; */
/*     for (int i = 1; i < argc; ++i) { */
/*         cout << "参数 " << i << ": " << argv[i] << endl; */
/*     } */
/*     return 0; */
/* } */
   
// ==================================== new 动态内存分配 ==========================================
/* #include <iostream> */
/* using namespace std; */

/* int main(int argc, char** argv) { */
/*     double *p = new double[10]; // 动态分配10个double的内存 */
/*     for (unsigned int i = 0; i < 10; i++) { */
/*         p[i] = (double)i; // 等价于 *(p+i) = i; */
/*     } */

/*     for (unsigned int i = 0; i < 10; i++) { */
/*         cout << p[i] << " "; */
/*     } */

/*     delete[] p; // 释放内存 */
/*     return 0; */
/* } */

// ======================================== 栈内存 =========================================
/* #include <iostream> */
/* using namespace std; */

/* int main() { */
/*     double p[10]; // 直接在栈上分配数组 */
/*     for (unsigned int i = 0; i < 10; i++) { */
/*         p[i] = (double)i; */
/*     } */

/*     for (unsigned int i = 0; i < 10; i++) { */
/*         cout << p[i] << " "; */
/*     } */
/*     return 0; */
/* } */

// =========================================== std::vector 自动管理 ========================================
#include <iostream>
#include <vector>
using namespace std;

int main() {
    vector<double> p(10); // 动态数组，自动管理内存
    for (unsigned int i = 0; i < 10; i++) {
        p[i] = (double)i;
    }

    // C++11 引入的 range-based for 循环（基于范围的 for 循环）
    for (auto val : p) { // 遍历每个元素
        cout << val << " "; // 输出该元素
    }
    cout << "\n";
    return 0;
}

















// ================================================== C++ class ================================================
/* #include <iostream> */
/* #include <fstream> */
/* using namespace std; */

/* class Vehicle { */
/*     public: */
/*         string brand = "Ford"; */
/*         void honk() { */
/*             cout << "Tuut, tuut\n"; */
/*         } */
/* }; */

/* class Car: public Vehicle { */
/*     public: */
/*         string model = "Mustang"; */
/* }; */


/* class myClass { */
/*     public: */
/*         void myFunction() { */
/*             cout << "Some content in parent class.\n"; */
/*         } */
/* }; */

/* class myOtherClass { */
/*     public: */
/*         void myOtherFunction() { */
/*             cout << "Some content in another class.\n"; */
/*         } */
/* }; */

/* class myChild: public myClass, public myOtherClass { */

/* }; */



/* class myGrandChild: public myChild { */

/* }; */


/* class Employee { */
/*     protected: */
/*         int salary; */
/* }; */

/* class Programmer: public Employee { */
/*     public: */
/*         int bonus; */
/*         void setSalary(int s) { */
/*             salary = s; */
/*         } */
/*         int getSalary() { */
/*             return salary; */
/*         } */
/* }; */


/* class Animal { */
/*     public: */
/*         void animalSound() { */
/*             cout << "the animal makes a sound\n"; */
/*         } */
/* }; */

/* class Pig: public Animal { */
/*     public: */
/*         void animalSound() { */
/*             cout << "the pig says: wee wee\n"; */
/*         } */
/* }; */

/* class Dog: public Animal { */
/*     public: */
/*         void animalSound() { */
/*             cout << "the dog says: bow bow\n"; */
/*         } */
/* }; */

/* class BaseParams { */
/*     public: */
/*         double getTime() const {return _t;}; */
/*         void setTime(double t0) {_t = t0;}; */

/*         double _t = 0; */
/* }; */
/* int main() { */
/*     /1* Car mycar; *1/ */
/*     /1* mycar.honk(); *1/ */
/*     /1* cout << mycar.brand + " " + mycar.model; *1/ */
/*     /1* myChild myobj; *1/ */
/*     /1* myobj.myFunction(); *1/ */
/*     /1* myobj.myOtherFunction(); *1/ */
/*     /1* Programmer guisong; *1/ */
/*     /1* guisong.setSalary(1000); *1/ */
/*     /1* guisong.bonus = 2000; *1/ */
/*     /1* cout << "salary:" << guisong.getSalary() << endl; *1/ */
/*     /1* cout << "salary:" << guisong.salary << endl; *1/ */
/*     /1* cout << "bonus:" << guisong.bonus << endl; *1/ */
/*     /1* Animal myAnimal; *1/ */
/*     /1* Dog myDog; *1/ */
/*     /1* Pig myPig; *1/ */
    
/*     /1* myAnimal.animalSound(); *1/ */
/*     /1* myDog.animalSound(); *1/ */
/*     /1* myPig.animalSound(); *1/ */
    
/*     /1* BaseParams myparam; *1/ */
/*     /1* myparam._t = 10; *1/ */
/*     /1* cout << myparam._t; *1/ */ 
/*     ofstream MyFile("filename.txt"); */
/*     MyFile << "Files can be tricky, but it is fun enough!"; */
/*     MyFile.close(); */
/*     return 0; */
/* } */

