FILE(REMOVE_RECURSE
  "CMakeFiles/mycode.out.dir/main.cc.o"
  "mycode.out.pdb"
  "mycode.out"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang CXX)
  INCLUDE(CMakeFiles/mycode.out.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
