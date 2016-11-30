FILE(REMOVE_RECURSE
  "CMakeFiles/test.out.dir/main.cc.o"
  "test.out.pdb"
  "test.out"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang CXX)
  INCLUDE(CMakeFiles/test.out.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
