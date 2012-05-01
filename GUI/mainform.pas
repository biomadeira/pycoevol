{*******************************************************************************
This file is part of the Pycoevol.
This work is public domain. Enjoy.
********************************************************************************
Author: Ludwig Krippahl
Date: 21.4.2012
Purpose:
Pycoevol GUI
Requirements:
Revisions:
To do:
*******************************************************************************}

unit mainform;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs, StdCtrls,
  Grids,Process, INIFiles;

type

  { TForm1 }

  TForm1 = class(TForm)
    StopPycoBt: TButton;
    PythonClEd: TEdit;
    Label10: TLabel;
    Label11: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    PycoMm: TMemo;
    OpenDialog1: TOpenDialog;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    File2BrowseBt: TButton;
    File2Ed: TEdit;
    Chain2Ed: TEdit;
    Id2Ed: TEdit;
    PsiblastCb: TComboBox;
    Label2: TLabel;
    AlignmentCb: TComboBox;
    CoevolutionCb: TComboBox;
    PycoFolderBrowseBt: TButton;
    ParamFileBrowseBt: TButton;
    File1BrowseBt: TButton;
    PycoFolderEd: TEdit;
    Label1: TLabel;
    ParamFileEd: TEdit;
    File1Ed: TEdit;
    Chain1Ed: TEdit;
    Id1Ed: TEdit;
    RunPyCoBt: TButton;
    SelectDirectoryDialog1: TSelectDirectoryDialog;
    procedure File1BrowseBtClick(Sender: TObject);
    procedure File2BrowseBtClick(Sender: TObject);
    procedure FormActivate(Sender: TObject);
    procedure FormClose(Sender: TObject; var CloseAction: TCloseAction);
    procedure FormCreate(Sender: TObject);
    procedure Label1Click(Sender: TObject);
    procedure ParamFileBrowseBtClick(Sender: TObject);
    procedure ParamFileEditClick(Sender: TObject);
    procedure PycoFolderBrowseBtClick(Sender: TObject);
    procedure RunPyCoBtClick(Sender: TObject);
    procedure StopPycoBtClick(Sender: TObject);
  private
    { private declarations }
    FTerminatePyCo:Boolean;
    Init:Boolean;
    function GetCommandLine:string;
    procedure LoadLists;
    procedure SaveConfiguration;
    procedure LoadConfiguration;
    procedure RunPycoevol;
  public
    { public declarations }
  end; 

var
  Form1: TForm1; 

implementation

{$R *.lfm}

{ TForm1 }

procedure TForm1.FormCreate(Sender: TObject);
begin
  LoadLists;
  Init:=True;
end;

procedure TForm1.Label1Click(Sender: TObject);
begin

end;

procedure TForm1.FormActivate(Sender: TObject);
begin
  if Init then
     begin
     LoadConfiguration;
     Init:=False;
     end;
end;

procedure TForm1.File2BrowseBtClick(Sender: TObject);
begin
  OpenDialog1.Filter:='PDB file|*.pdb|Sequence|*.fasta|Any|*.*';
  if OpenDialog1.Execute then
     File2Ed.Text:=OpenDialog1.FileName;
end;

procedure TForm1.File1BrowseBtClick(Sender: TObject);
begin
  OpenDialog1.Filter:='PDB file|*.pdb|Sequence|*.fasta|Any|*.*';
  if OpenDialog1.Execute then
     File1Ed.Text:=OpenDialog1.FileName;
end;

procedure TForm1.FormClose(Sender: TObject; var CloseAction: TCloseAction);
begin
  SaveConfiguration;
end;

procedure TForm1.ParamFileBrowseBtClick(Sender: TObject);
begin
  OpenDialog1.Filter:='Parameter file|*.config|Any|*.*';
  if OpenDialog1.Execute then
     ParamFileEd.Text:=OpenDialog1.FileName;
end;

procedure TForm1.ParamFileEditClick(Sender: TObject);

var proc:TProcess;

begin
  proc:=TProcess.Create(nil);
  proc.CommandLine:='"'+ParamFileEd.Text+'"';
  proc.Execute;
  proc.Free;
end;

procedure TForm1.PycoFolderBrowseBtClick(Sender: TObject);
begin
  if SelectDirectoryDialog1.Execute then
    PycoFolderEd.Text:=SelectDirectoryDialog1.FileName;
end;

procedure TForm1.RunPyCoBtClick(Sender: TObject);

var oc:string;

begin
  RunPyCoBt.Enabled:=False;
  StopPyCoBt.Enabled:=True;
  FTerminatePyco:=False;
  oc:=RunPyCoBt.Caption;
  RunPyCoBt.Caption:='Busy...';
  Application.ProcessMessages;
  try
    RunPycoevol;
  finally
    RunPyCoBt.Caption:=oc;
    RunPyCoBt.Enabled:=True;
  end;
end;

procedure TForm1.StopPycoBtClick(Sender: TObject);
begin
  StopPyCoBt.Enabled:=False;
  FTerminatePyco:=True;
end;

function TForm1.GetCommandLine: string;
begin

  //TODO: Check for spaces in parameters??

  Result:=PythonClEd.Text+' "'+PycoFolderEd.Text+PathDelim+'Pycoevol.py"';
  if (File1Ed.Text<>'') and (File2Ed.Text<>'') then
    Result:=Result+' "'+File1Ed.Text+'" "'+File2Ed.Text+'"';
  if (Chain1Ed.Text<>'') and (Chain2Ed.Text<>'') then
    Result:=Result+' -x'+Chain1Ed.Text+' -x'+Chain2Ed.Text;
  if (Id1Ed.Text<>'') and (Id2Ed.Text<>'') then
    Result:=Result+' -i'+Id1Ed.Text+' -i'+Id2Ed.Text;
  Result:=Result+' -b'+PsiblastCb.Text+' -a'+AlignmentCb.Text+' -c'+CoevolutionCb.Text;
  if (ParamFileEd.text<>'') then
  Result:=Result+' -p"'+ParamFileEd.text+'"';
end;

procedure TForm1.LoadLists;

var
  sl:TStringList;
  f:Integer;
  s:string;
  currbox:TComboBox;

begin
  sl:=TStringList.Create;
  sl.LoadFromFile('optionlists.txt');
  for f:=0 to sl.Count-1 do
    begin
    s:=sl.Strings[f];
    if s='**coevolution' then currbox:=CoevolutionCB
    else if s='**alignment' then currbox:=AlignmentCB
    else if s='**psiblast' then currbox:=PsiblastCB
    else
      currbox.Items.Add(s);
    end;
  CoevolutionCb.ItemIndex:=0;
  AlignmentCb.ItemIndex:=0;
  PsiblastCb.ItemIndex:=0;

  sl.Free;
end;

procedure TForm1.SaveConfiguration;

var
  cfg:string;
  ini:TIniFile;

begin
  cfg := GetAppConfigFile(False);
  if not DirectoryExists(ExtractFileDir(cfg)) then
      CreateDir(ExtractFileDir(cfg));
  ini:=TiniFile.Create(cfg);
  ini.WriteString('Form','PycoevolFolder',PycoFolderEd.Text);
  ini.WriteString('Form','ParametersFile',ParamFileEd.Text);
  ini.WriteString('Form','File1',File1Ed.Text);
  ini.WriteString('Form','File2',File2Ed.Text);
  ini.WriteString('Form','Psiblast',PsiblastCb.Text);
  ini.WriteString('Form','Alignment',AlignmentCB.Text);
  ini.WriteString('Form','Coevolution',CoevolutionCB.Text);
  ini.WriteString('Form','Python',PythonClEd.Text);
  ini.UpdateFile;
  ini.Free;
end;

procedure TForm1.LoadConfiguration;

var
   cfg:string;
   ini:TIniFile;

begin
  cfg := GetAppConfigFile(False);
  if FileExists(cfg) then
    begin
    ini:=TIniFile.Create(cfg);
    PycoFolderEd.Text:=ini.ReadString('Form','PycoevolFolder','');
    ParamFileEd.Text:=ini.ReadString('Form','ParametersFile','');
    File1Ed.Text:=ini.ReadString('Form','File1','');
    File2Ed.Text:=ini.ReadString('Form','File2','');
    PythonClEd.Text:=ini.ReadString('Form','Python','');

    PsiblastCb.ItemIndex:=PsiblastCb.Items.IndexOf(ini.ReadString('Form','Psiblast',''));
    if PsiBlastCB.ItemIndex<0 then PsiBlastCb.ItemIndex:=0;
    AlignmentCB.ItemIndex:=AlignmentCB.Items.IndexOf(ini.ReadString('Form','Alignment',''));
    if AlignmentCB.ItemIndex<0 then AlignmentCB.ItemIndex:=0;
    CoevolutionCB.ItemIndex:=CoevolutionCB.Items.IndexOf(ini.ReadString('Form','Coevolution',''));
    if CoevolutionCB.ItemIndex<0 then CoevolutionCB.ItemIndex:=0;

    ini.Free;
    end;
end;

procedure TForm1.RunPycoevol;


var
  ebytes,nbytes: LongInt;
  proc:TProcess;
  cl:string;

procedure RefreshOutput;

var
  s:string;

begin
  nbytes := proc.Output.NumBytesAvailable;
  while nbytes > 0 do
     begin
     SetLength(s,nbytes);
     proc.Output.Read(s[1], nbytes);
     nbytes := proc.Output.NumBytesAvailable;
     PycoMm.Lines.Add(s);
     end;
  ebytes := proc.Stderr.NumBytesAvailable;
  if ebytes>0 then PycoMm.Lines.Add('*** ERROR ***');
  while ebytes > 0 do
     begin
     SetLength(s,ebytes);
     proc.Stderr.Read(s[1], ebytes);
     ebytes := proc.Stderr.NumBytesAvailable;
     PycoMm.Lines.Add(s);
     end;
   PycoMm.SelStart := Length(PycoMm.Lines.Text)-1;
   PycoMm.SelLength:=0;
   PycoMm.SetFocus;
   Application.ProcessMessages;
   if FTerminatePyCo then proc.Active:=False;
end;

begin
  SetCurrentDir(PyCoFolderEd.Text);
  cl:=GetCommandLine;
  proc:=TProcess.Create(nil);
  proc.CommandLine := cl;
  proc.Options := [poUsePipes];
  PycoMm.Lines.Add(cl);
  Application.ProcessMessages;
  proc.Execute;
  while proc.Running do
  begin
    Application.ProcessMessages;
    RefreshOutput;
    Sleep(500);
  end;
  RefreshOutput;
  proc.Free;
  if FTerminatePyCo then PycoMm.Lines.Add('Terminated by user')
  else PycoMm.Lines.Add('Done');
  PycoMm.SelStart := Length(PycoMm.Lines.Text)-1;
  PycoMm.SelLength:=0;
  PycoMm.SetFocus;

end;

end.

