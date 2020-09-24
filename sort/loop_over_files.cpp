#include "TSystem.h"

void loop_over_files()
{
    const char* ext = ".root";
    const char* inDir = "/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C";
    char* dir = gSystem->ExpandPathName(inDir);
    void* dirp = gSystem->OpenDirectory(dir);
    const char* entry;
    const char* filename[200];
    Int_t n = 0;
    TString str;
    while((entry = (char*)gSystem->GetDirEntry(dirp))) {
        str = entry;
        if(str.EndsWith(ext))
            filename[n++] = gSystem->ConcatFileName(dir, entry);
    }
    for (Int_t i = 0; i < n; i++){
        Printf("file -> %s", filename[i]);
    }
    gSystem->FreeDirectory(dirp);
}
