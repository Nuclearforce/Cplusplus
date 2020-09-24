#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <iostream>
#include <vector>
#include <TRandom3.h>
#include <cstdio>

//start 'global'
TRandom3 *r = new TRandom3();
TFile f_hist("../histogram/Gamma_calibration/fixed_Eu_calibrated_after_background_subtracted_exp_6March.root", "recreate");
//TFile f_hist("../histogram/13C_uncalibrated_before_experiment.root", "recreate");
TChain DataChain("data");

TTreeReader myReader(&DataChain);
TTreeReaderArray<unsigned int> TRenergy(myReader, "energy");
TTreeReaderArray<unsigned int> TRtime(myReader, "time");
TTreeReaderArray<unsigned int> TRmodchan(myReader, "modchan");

TChain Background("data");
TTreeReader myReader_background(&Background);
TTreeReaderArray<unsigned int> background_energy(myReader_background, "energy");
TTreeReaderArray<unsigned int> background_time(myReader_background, "time");
TTreeReaderArray<unsigned int> background_modchan(myReader_background, "modchan");
//histograms
TH1F* hPixieADC[4][16];
TH1F* hPixieADC_subtracted[4][16];
TH1F* backgroundADC[4][16];
TH1F* Clover3;
TH1F* Clover1;
//calibration
Double_t calib_a[4][16];
Double_t calib_b[4][16];
Double_t calib_c[4][16];
Double_t calib_d[4][16];
Double_t calib_e[4][16];

//Doppler correction
Double_t beta = 0.066;
Double_t theta[4][16];
Double_t phi[4][16];

//timing
Double_t sector_time;
Double_t ring_time;

//function declaration
void second_loop();
void init_histos();
void read_calibrate_file();
void make_chain();
void first_loop();
Double_t do_calibrate(Double_t raw_energy, Int_t module, Int_t channel);
//Double_t do_doppler(Double_t en,Int_t gamma_mod, Int_t gamma_chn);
//Double_t do_beta(Double_t ring_en);

void three_leaves_gamma_calib()
{
    
    make_chain();//add trees to chain
    init_histos(); //initialise histograms
    read_calibrate_file(); //read calibration file
    
    // start timer for progress bar
    std::clock_t start;
    Double_t duration;
    Long64_t timer = 0;
    //Long64_t max_event = 1e6;
    Long64_t max_event = DataChain.GetEntries();
    Long64_t background_max_event = Background.GetEntries();
    start = std::clock();
    
    // Loop over all entries of the TChain.
    while (myReader.Next()) {
        first_loop(); //fills singles histograms for now
        if(timer%4000000==0){
            cout << "Event #: " << timer  << " Progress: " << 100.0*timer/(max_event+background_max_event) << " [%] \r";
            cout.flush();
        }
        timer++;
    }
    cout<<"spam1"<<endl;
    //Loop over background
    while (myReader_background.Next()) {
        
        second_loop(); //fills singles histograms for now
        if(timer%4000000==0){
        cout << "Event #: " << timer  << " Progress: " << 100.0*timer/(max_event+background_max_event) << " [%] \r";
        cout.flush();
        }
        timer++;
    }
    cout<<"spam2"<<endl;
    Int_t j=0;
    Int_t i=0;
    for(j=0;j<4;j++){
        
        for(i=0;i<16;i++){
            hPixieADC_subtracted[j][i] -> Add(hPixieADC[j][i],1);
            hPixieADC_subtracted[j][i] -> Add(backgroundADC[j][i],-0.350194552);
        }
    }
    
    Clover1->Add(hPixieADC[1][0],1);
    Clover1->Add(hPixieADC[1][1],1);
    Clover1->Add(hPixieADC[1][2],1);
    Clover1->Add(hPixieADC[1][3],1);
    Clover3->Add(hPixieADC[3][0],1);
    Clover3->Add(hPixieADC[3][1],1);
    Clover3->Add(hPixieADC[3][2],1);
    Clover3->Add(hPixieADC[3][3],1);
    cout<<"spam3"<<endl;
    f_hist.Write();
    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    cout << "Timer: " << duration << " [s] to analyze " << max_event << " events." << endl;
    f_hist.Close();
    
}

void second_loop(){
    Int_t entry_id_back;
    Int_t entry_size_back;
    
    Double_t time_back;
    Double_t en_back;
    Double_t calib_en_back;
    
    for (entry_id_back = 0, entry_size_back = background_energy.GetSize(); entry_id_back < entry_size_back; ++entry_id_back) {
        
        Int_t modchan=background_modchan[entry_id_back];
        Int_t mod=modchan/16;
        Int_t chn=modchan%16;
        time_back=background_time[entry_id_back];
        en_back=background_energy[entry_id_back];
        calib_en_back=do_calibrate(en_back,mod,chn);
        backgroundADC[mod][chn]->Fill(calib_en_back);
        
    }

}

void first_loop(){
    
    Int_t entry_id;
    Int_t entry_size;
    Int_t sector_count = 0;
    Int_t ring_count = 0;
    Double_t time;
    
    Double_t en;
    Double_t calib_en;
    
    for (entry_id = 0, entry_size = TRenergy.GetSize(); entry_id < entry_size; ++entry_id) {
        
        Int_t modchan=TRmodchan[entry_id];
        Int_t mod=modchan/16;
        Int_t chn=modchan%16;
        time=TRtime[entry_id];
        en=TRenergy[entry_id];
        calib_en=do_calibrate(en,mod,chn);
        hPixieADC[mod][chn]->Fill(calib_en);
        //if ((mod==3)&&(chn==0)){
        //cout<<entry_size<< " " << entry_id << " " <<mod<< " " << chn << " " << time << " " << en <<endl;
        //}
    }
    
    
}


void read_calibrate_file() {
    Double_t a = 0.0;
    Double_t b = 1.0;
    Double_t c = 1.0;
    Double_t d = 0.0;
    Double_t e = 0.0;
    Double_t f = 0.0;
    Double_t g = 0.0;
    FILE *calib = fopen("Calibration_after_exp.txt","r");
        for(int m=0;m<4;m++){
        for(int n=0;n<16;n++){ //The format is: offset(a) slope^1(b) slope(c)^2 slope(f)^3 slope(g)^4 phi(d) theta(e)
            fscanf(calib, "%lf %lf %lf %lf %lf %lf %lf", &a, &b, &c, &f, &g, &d, &e);
            calib_a[m][n]=a;
            calib_b[m][n]=b;
            calib_c[m][n]=c;
            calib_d[m][n]=f;
            calib_e[m][n]=g;
            phi[m][n]=d;
            theta[m][n]=e;
        }
    }
    
}

Double_t do_calibrate(Double_t raw_energy, Int_t module, Int_t channel){
    Double_t random_num;
    Double_t cal_energy;
    Double_t offset = calib_a[module][channel];
    Double_t slope1 = calib_b[module][channel];
    Double_t slope2 = calib_c[module][channel];
    Double_t slope3 = calib_d[module][channel];
    Double_t slope4 = calib_e[module][channel];
    random_num=r->Uniform(-0.5,0.5);
    cal_energy=(raw_energy+random_num)*(raw_energy+random_num)*(raw_energy+random_num)*(raw_energy+random_num)*slope4+(raw_energy+random_num)*(raw_energy+random_num)*(raw_energy+random_num)*slope3+(raw_energy+random_num)*(raw_energy+random_num)*slope2 + (raw_energy+random_num)*slope1 + offset;
   /*if (module==2&&channel==4&&raw_energy>8300&&raw_energy<8500){
        cout << raw_energy << " " << slope1 << " " << slope2 << " " << slope3 << " " << slope4 << " "  << cal_energy <<endl;
    }*/
    
    return cal_energy;
}

/*
 Double_t do_beta(Double_t ring_en){
 Double_t beta;
 beta=sqrt(((ring_en/1000.0)*2.0)/(931.5*14.0));
 return beta;
 }
 
 Double_t do_doppler(Double_t en, Int_t gamma_mod, Int_t gamma_chn){
 Double_t dop_cor_en;
 Double_t cos_alpha;
 cos_alpha=sin(theta[ring_mod][ring_chn])*sin(theta[gamma_mod][gamma_chn])*cos(phi[sector_mod][sector_chn]-phi[gamma_mod][gamma_chn]) + cos(theta[ring_mod][ring_chn])*cos(theta[gamma_mod][gamma_chn]);
 dop_cor_en=(en)*(1.0-beta*cos_alpha)/(sqrt(1-beta*beta));
 return dop_cor_en;
 }
 */

void init_histos() {
    //calibrated
    /*Int_t gamma_bins=10000;
     Int_t gamma_min=0;
     Int_t gamma_max=9999;
     Int_t particle_bins=36000;
     Int_t particle_min=0;
     Int_t particle_max=35999;*/
    
    //uncalibrated
    Int_t gamma_bins=32768;
    Int_t gamma_min=0;
    Int_t gamma_max=32768;
    Int_t particle_bins=53000;
    Int_t particle_min=0;
    Int_t particle_max=53000;
    
    
    /*
     hPixieADC[1][0] = new TH1F("CL1_Red","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[1][1] = new TH1F("CL1_Green","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[1][2] = new TH1F("CL1_Blue","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[1][3] = new TH1F("CL1_Black","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[3][0] = new TH1F("CL3_Red","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[3][1] = new TH1F("CL3_Green","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[3][2] = new TH1F("CL3_Blue","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[3][3] = new TH1F("CL3_Black","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[2][0] = new TH1F("LaBr3_1","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[2][1] = new TH1F("LaBr3_2","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[2][2] = new TH1F("LaBr3_3","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[2][3] = new TH1F("LaBr3_4","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[2][4] = new TH1F("LaBr3_5","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[2][5] = new TH1F("LaBr3_6","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[3][5] = new TH1F("FSU7","",gamma_bins,gamma_min,gamma_max);
     hPixieADC[0][0] = new TH1F("Sector32_1","",particle_bins,particle_min,particle_max);
     hPixieADC[0][1] = new TH1F("Sector2_3","",particle_bins,particle_min,particle_max);
     hPixieADC[0][2] = new TH1F("Sector4_5","",particle_bins,particle_min,particle_max);
     hPixieADC[0][3] = new TH1F("Sector6_7","",particle_bins,particle_min,particle_max);
     hPixieADC[0][4] = new TH1F("Sector8_9","",particle_bins,particle_min,particle_max);
     hPixieADC[0][5] = new TH1F("Sector10_11","",particle_bins,particle_min,particle_max);
     hPixieADC[0][6] = new TH1F("Sector12_13","",particle_bins,particle_min,particle_max);
     hPixieADC[0][7] = new TH1F("Sector14_15","",particle_bins,particle_min,particle_max);
     hPixieADC[0][8] = new TH1F("Sector16_17","",particle_bins,particle_min,particle_max);
     hPixieADC[0][9] = new TH1F("Sector18_19","",particle_bins,particle_min,particle_max);
     hPixieADC[0][10] = new TH1F("Sector20_21","",particle_bins,particle_min,particle_max);
     hPixieADC[0][11] = new TH1F("Sector22_23","",particle_bins,particle_min,particle_max);
     hPixieADC[0][12] = new TH1F("Sector24_25","",particle_bins,particle_min,particle_max);
     hPixieADC[0][13] = new TH1F("Sector26_27","",particle_bins,particle_min,particle_max);
     hPixieADC[0][14] = new TH1F("Sector28_29","",particle_bins,particle_min,particle_max);
     hPixieADC[0][15] = new TH1F("Sector30_31","",particle_bins,particle_min,particle_max);
     hPixieADC[1][8] = new TH1F("R1","",particle_bins,particle_min,particle_max);
     hPixieADC[1][9] = new TH1F("R2","",particle_bins,particle_min,particle_max);
     hPixieADC[1][10] = new TH1F("R3","",particle_bins,particle_min,particle_max);
     hPixieADC[1][11] = new TH1F("R4","",particle_bins,particle_min,particle_max);
     hPixieADC[1][12] = new TH1F("R5","",particle_bins,particle_min,particle_max);
     hPixieADC[1][13] = new TH1F("R6","",particle_bins,particle_min,particle_max);
     hPixieADC[1][14] = new TH1F("R7","",particle_bins,particle_min,particle_max);
     hPixieADC[1][15] = new TH1F("R8","",particle_bins,particle_min,particle_max);
     hPixieADC[2][8] = new TH1F("R9","",particle_bins,particle_min,particle_max);
     hPixieADC[2][9] = new TH1F("R10","",particle_bins,particle_min,particle_max);
     hPixieADC[2][10] = new TH1F("R11","",particle_bins,particle_min,particle_max);
     hPixieADC[2][11] = new TH1F("R12","",particle_bins,particle_min,particle_max);
     hPixieADC[2][12] = new TH1F("R13","",particle_bins,particle_min,particle_max);
     hPixieADC[2][13] = new TH1F("R14","",particle_bins,particle_min,particle_max);
     hPixieADC[2][14] = new TH1F("R15","",particle_bins,particle_min,particle_max);
     hPixieADC[2][15] = new TH1F("R16","",particle_bins,particle_min,particle_max);
     hPixieADC[3][8] = new TH1F("R17","",particle_bins,particle_min,particle_max);
     hPixieADC[3][9] = new TH1F("R18","",particle_bins,particle_min,particle_max);
     hPixieADC[3][10] = new TH1F("R19","",particle_bins,particle_min,particle_max);
     hPixieADC[3][11] = new TH1F("R20","",particle_bins,particle_min,particle_max);
     hPixieADC[3][12] = new TH1F("R21","",particle_bins,particle_min,particle_max);
     hPixieADC[3][13] = new TH1F("R22","",particle_bins,particle_min,particle_max);
     hPixieADC[3][14] = new TH1F("R23","",particle_bins,particle_min,particle_max);
     hPixieADC[3][15] = new TH1F("R24","",particle_bins,particle_min,particle_max);
     hPixieADC[1][5] = new TH1F("empty1","",2,0,1);
     hPixieADC[1][6] = new TH1F("empty2","",2,0,1);
     hPixieADC[1][7] = new TH1F("empty3","",2,0,1);
     hPixieADC[2][6] = new TH1F("empty4","",2,0,1);
     hPixieADC[2][7] = new TH1F("empty5","",2,0,1);
     hPixieADC[3][7] = new TH1F("empty6","",2,0,1);
     hPixieADC[3][6] = new TH1F("FSU7_BGO","",2,0,1);
     hPixieADC[1][4] = new TH1F("CL1_BGO","",2,0,1);
     hPixieADC[3][4] = new TH1F("CL3_BGO","",2,0,1);
     */
    Clover1 = new TH1F("Clover_1","",particle_bins,particle_min,particle_max);
    Clover3 = new TH1F("Clover_3","",particle_bins,particle_min,particle_max);
    
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            hPixieADC[j][i] = new TH1F(Form("ADC%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            hPixieADC_subtracted[j][i] = new TH1F(Form("subtracted_ADC%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            backgroundADC[j][i] = new TH1F(Form("ADC_background%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    
    
}

void make_chain(){
    
    //////14C/////////////
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0678-00.root");
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0806-00.root");
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0915-00.root");
    
    
    ///////13C before experiment////////////
    /*DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0651-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0653-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0653-01.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0654-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0654-01.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0655-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0656-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0657-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0657-01.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0658-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0658-01.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0659-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0660-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0660-01.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0661-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0661-01.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0662-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0662-01.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f31d2-44d7-8418-e1a6319e7468/13C/run-0663-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0664-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0664-01.root");*/
    
    ///////////Si calibration/////////
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0819-00.root");//38 MeV
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0820-00.root");//30 MeV
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0821-00.root");//52.6MeV
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0817-00.root");//45 MeV
    
    ///////////13C after experiment/////////
    /*DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0918-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0918-01.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0918-02.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0918-03.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0918-04.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0919-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0919-01.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0919-02.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0919-03.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0920-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0920-01.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0920-02.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0921-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0921-01.root");*/
    
    /////////////152Eu//////////////////
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0668-00.root");//58min run
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0926-00.root"); //3 hours
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0926-01.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0926-02.root");
    ///////////background runs before//////////////
    /*Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0670-00.root");//8.48hours
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0670-01.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0671-00.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0672-00.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0672-01.root");*/
    ///////////background runs after//////////////
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-00.root");//08:33:38
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-01.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-02.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-03.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-04.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-05.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-06.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-07.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-08.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-09.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-10.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-11.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-12.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-13.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-14.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-15.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-16.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-17.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-18.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-19.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-20.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-21.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-22.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-23.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-24.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-25.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-26.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-27.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-28.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-29.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-30.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-31.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-32.root");
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-33.root");

    //////////////Ga66 after experiment beamoff/////////
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0923-00.root");//10min
    /*DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0924-00.root");//1 day and 33min
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0924-01.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0924-02.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0924-03.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0924-04.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0924-05.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0924-06.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0924-07.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0924-08.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0924-09.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0924-10.root");*/

    
    //////////////Ga66 after experiment beam on/////////
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0922-00.root");
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0922-01.root");
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0922-02.root");
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0922-03.root");
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0922-04.root");
    
    //////////////Ga66 before experiment beamoff/////////
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0557.root");//18:13min
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0558.root");//17:45 hours
    
    //////////////Ga66 before experiment beam on/////////
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0556.root");
    
    
}
