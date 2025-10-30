# My-projects

#for this to work following is required 
#compiler: mnvs 19.4xx  
#Cmake 3.18 or higher 
#SFML 2.6.1 or 2.6.2 
#dowload imgui and imgui-sfml from GitHub 
_________________________________________________________________
#don't put readme file here

#file structure
BlackHoleProject/
├── CMakeLists.txt              
├── main.cpp                    
├── stb_image_write.h          

├── imgui/                      # <--- Put ImGui files here.
│   ├── imgui.h                 
│   ├── imgui.cpp
│   ├── imgui_draw.cpp
│   ├── imgui_widgets.cpp
│   └── imgui_internal.h        
│   └── imconfig.h              
│   └── ImFontAtlas.cpp         
│   └── ImFontAtlas.h
│   └── ImSequencer.h           
|    └── ImVec2.h               # (and so on... download ALL imgui source files)         
|
└── imgui-sfml/                 # <--- Create this FOLDER.
    ├── imgui-SFML.h            # <--- Put ImGui-SFML files here.
    ├── imgui-SFML.cpp
    └── imgui-SFML_export.h     
    └── ImGuiExt.h              # (and so on... download ALL imgui-sfml source files)

____________________________________________________________________________________

#open project file on terminal then... 
#1.mkdir build   (makes a build file)
#2. $env:SFML_DIR=" path_of_smfl"  
{#example: $env:SFML_DIR="C:\Users\swastik\Downloads\SFML-3.0.2\lib\cmake\SFML" }[this sets envormental variable ]  
#3.cd build  
#4.cmake ..   
#5.cmake --build .    
#6../bin/"paste the .exe file adress here found in bin folder"    this is to run project

#warning!!!: mismatch in compiler version, cmake version or sfml version can cause error regarding imgui syntax
