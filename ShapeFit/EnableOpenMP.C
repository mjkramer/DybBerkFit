TString arrstr(TObjArray& array, size_t i)
{
  return ((TObjString*)(array[i]))->GetString();
}

void EnableOpenMP()
{
  TString cmd = gSystem->GetMakeSharedLib();
  TObjArray& parts = *cmd.Tokenize(";");
  TString cd = arrstr(parts, 0);
  TString compile = arrstr(parts, 1);
  TString link = arrstr(parts, 2);
  const char* newCmd = Form("%s ; %s -fopenmp ; %s -fopenmp",
                            cd.Data(), compile.Data(), link.Data());
  gSystem->SetMakeSharedLib(newCmd);
}
