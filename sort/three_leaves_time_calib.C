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
TFile f_hist("../histogram/time_calibration_check_week_three.root", "recreate");
//TFile f_hist("../histogram/13C_uncalibrated_before_experiment.root", "recreate");
TChain DataChain("data");
TTreeReader myReader(&DataChain);
TTreeReaderArray<unsigned int> TRenergy(myReader, "energy");
TTreeReaderArray<unsigned int> TRtime(myReader, "time");
TTreeReaderArray<unsigned int> TRmodchan(myReader, "modchan");
//histograms
TH1F* hPixieADC[4][16];
TH1F* LaBr_particle_gamma_time;
TH1F* Ge_particle_gamma_time;
TH1F* S3_ring_time;
TH1F* S3_sector_time;
TH1F* S3_ring_sector_time;
TH2F* Gamma_time;
TH2F* silicon_time;

//calibration
Double_t calib_a[4][16];
Double_t calib_b[4][16];

//Doppler correction
Double_t beta = 0.066;
Double_t theta[4][16];
Double_t phi[4][16];

//timing
Double_t sector_time;
Double_t ring_time;

//function declaration
void init_histos();
void read_calibrate_file();
void make_chain();
void first_loop();
void second_loop();
Double_t do_calibrate(Double_t raw_energy, Int_t module, Int_t channel);
//Double_t do_doppler(Double_t en,Int_t gamma_mod, Int_t gamma_chn);
//Double_t do_beta(Double_t ring_en);

void three_leaves()
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
    start = std::clock();
   
    // Loop over all entries of the TChain.
    while (myReader.Next()) {
        
        first_loop(); //fills singles histograms for now
        second_loop(); //fills histograms that need gates
        if(timer%1000000==0){
            cout << "Event #: " << timer  << " Progress: " << 100.0*timer/max_event << " [%] \r";
            cout.flush();
        }
        timer++;
    }
    
    f_hist.Write();
    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    cout << "Timer: " << duration << " [s] to analyze " << max_event << " events." << endl;
    f_hist.Close();
    
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
        
        
        if ((mod==0)&&(sector_count==0)){
            S3_sector_time -> Fill(time+200);
            sector_time=time;
            sector_count++;
            
        }
        
        if ( (((mod==1)&&(chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15))||(((mod==2)&&(chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15))) || (((mod==3)&&(chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)))) && (ring_count==0)){
            
            
            S3_ring_time -> Fill(time+200);
            sector_time=time;
            ring_count++;
            
        }
        //if ((mod==3)&&(chn==0)){
        //cout<<entry_size<< " " << entry_id << " " <<mod<< " " << chn << " " << time << " " << en <<endl;
        //}
    }
    
    
}

void second_loop(){
    
    Int_t entry_id;
    Int_t entry_size;
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
        S3_ring_sector_time->Fill(sector_time-ring_time+200);
        
        if (mod==0){
            silicon_time->Fill(time+200,chn);
        }
        if(mod==1){
            if (chn<8){
                
                Gamma_time->Fill(sector_time-time+200,chn);
                
            }
            else{
                
                silicon_time->Fill(time+200,chn+8);
                
            }
        }
        if(mod==2){
            if (chn<8){
                
                Gamma_time->Fill(sector_time-time+200,chn+8);
                
            }
            else{
                
                silicon_time->Fill(time+200,chn+16);
                
            }
        }
        if(mod==3){
            if (chn<8){
                
                Gamma_time->Fill(sector_time-time+200,chn+16);
                
            }
            else{
                
                silicon_time->Fill(time+200,chn+24);
                
            }
        }

        if( (mod==2)&&(chn==0||chn==1||chn==2||chn==3||chn==4||chn==5) ){
            
            LaBr_particle_gamma_time->Fill(sector_time-time+200);
            
        }
        
        
        if(((mod==1)&&(chn==0||chn==1||chn==2||chn==3))||((mod==3)&&(chn==0||chn==1||chn==2||chn==3))){
            
            Ge_particle_gamma_time->Fill(sector_time-time+200);
            
        }

        
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
    FILE *calib = fopen("Calibration_empty.txt","r");
    //FILE *calib = fopen("Calibration_14C.txt","r");
    for(int m=0;m<4;m++){
        for(int n=0;n<16;n++){ //The format is: offset(a) slope(b) phi(c) theta(d)
            fscanf(calib, "%lf %lf %lf %lf %lf", &a, &b, &c, &d, &e);
            calib_a[m][n]=a;
            calib_b[m][n]=b;
            calib_c[m][n]=c;
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
    random_num=r->Uniform(-0.5,0.5);
    cal_energy=(raw_energy+random_num)*slope2*slope2 + (raw_energy+random_num)*slope1 + offset;
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
    Int_t particle_bins=32768;
    Int_t particle_min=0;
    Int_t particle_max=32768;
    
    LaBr_particle_gamma_time = new TH1F("LaBr_particle_gamma_time","",400,0,399);
    Ge_particle_gamma_time = new TH1F("Ge_particle_gamma_time","",400,0,399);
    S3_ring_time = new TH1F("S3_ring_time","",400,0,399);
    S3_sector_time = new TH1F("S3_sector_time","",400,0,399);
    S3_ring_sector_time = new TH1F("S3_ring_sector_time","",400,0,399);
    Gamma_time = new TH2F("Gamma_time","",20,0,399,24,0,23);
    silicon_time = new TH2F("silicon_time","",20,0,399,32,0,31);
    
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
    Int_t j=0;
    Int_t i=0;
    for(j=0;j<4;j++){
        
        for(i=0;i<16;i++){
            hPixieADC[j][i] = new TH1F(Form("ADC%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
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
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0817-00.root");//45 MeV
    
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
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0668-00.root");//58min run
    
    ///////////background runs//////////////
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0670-00.root");//2h44min
    
    //////////////Ga66 before experiment beamoff/////////
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0557.root");
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0558.root");
    
    //////////////Ga66 before experiment beam on/////////
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0556.root");

    
}
