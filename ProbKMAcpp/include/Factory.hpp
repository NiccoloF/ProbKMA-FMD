#ifndef __BUILDER_HPP__
#define __BUILDER_HPP__

#include <unordered_map>
#include <string>
#include <functional>
#include <memory>

// Factory class
namespace util{
 template<typename B>
 class SharedFactory
 {
   public :
     
     using registry_map = std::unordered_map<std::string,std::function<std::unique_ptr<B>()>>;
     
     registry_map map;
     
     // use this to instantiate the proper Derived class
     std::unique_ptr<B> instantiate (const std::string& name)
       {
         auto it = map.find(name);
         return it == map.end () ? nullptr : (it->second)();
         }
    
    template <typename D, typename... Args>
    void FactoryRegister(const std::string& name, Args&&... args) {
      map[name] = [args = std::forward_as_tuple(std::forward<Args>(args)...)]() mutable 
        {
          return std::apply([](auto&&... a) 
          {return std::make_unique<D>(std::forward<decltype(a)>(a)...);}, std::move(args));
        };
    
     };
 };

}


#endif //__BUILDER_HPP__