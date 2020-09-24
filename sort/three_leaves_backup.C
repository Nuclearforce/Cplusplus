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
TFile f_hist("../histogram/efficiency/after/13C_doppler_calib_test.root", "recreate");
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

//extra
Bool_t for_testing = true;

//gate
Bool_t gate_on_excitation_spectra;
Bool_t particle_gamma;
Bool_t energy_sharing_ring;
Bool_t energy_sharing_sector;
Bool_t ring_sector_time_gate;
Double_t ring_time;
Double_t sector_time;
Double_t gamma_time;
Bool_t bad_gamma;

//lower and upper limit for excitation spectrum gate
Double_t lower_lim = 1600;
Double_t upper_lim = 1683;

//histograms
TH1F* hPixieADC[4][16];
TH1F* hPixieADC_ab[4][16];
TH1F* hPixieADC_subtracted[4][16];
TH1F* hPixieADC_pg[4][16];
TH1F* backgroundADC[4][16];
TH1F* time_Ge;
TH1F* time_La;
TH1F* time_La_army;
TH1F* pg_time;
TH1F* ring_sector_time;
TH1F* time_background[4][16];
TH1F* time_Si;

//calibration
Double_t calib_a[4][16];
Double_t calib_b[4][16];
Double_t calib_c[4][16];

//Doppler correction
Double_t beta = 0.017947;
Double_t theta[4][16];
Double_t phi[4][16];
Int_t ring_mod;
Int_t ring_chn;
Int_t sector_mod;
Int_t sector_chn;

//energy_sharing
Double_t energy_sharing_ring_en=0;
Double_t energy_sharing_ring_chn_use=0;
Double_t energy_sharing_ring_mod_use=0;
Double_t energy_sharing_ring_chn_ignore=0;
Double_t energy_sharing_ring_mod_ignore=0;
Double_t energy_sharing_sector_en=0;
Double_t energy_sharing_sector_chn_use=0;
Double_t energy_sharing_sector_chn_ignore=0;
Bool_t bad_particle;

//Add-back
Double_t ab_use_energy_cl1;
Double_t ab_use_crystal_cl1;
Double_t ab_ignore_crystal_cl1;
Bool_t ab_event_cl1=false;
Double_t ab_use_energy_cl3;
Double_t ab_use_crystal_cl3;
Double_t ab_ignore_crystal_cl3;
Bool_t ab_event_cl3=false;
Double_t test_counter=0;

//function declaration
void background_loop();
void init_histos();
void read_calibrate_file();
void make_chain();
void first_loop();
void second_loop();
Double_t do_calibrate(Double_t raw_energy, Int_t module, Int_t channel);
Double_t do_doppler(Double_t en,Int_t gamma_mod, Int_t gamma_chn);
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
    Long64_t max_event_data = DataChain.GetEntries();
    Long64_t max_event_background = Background.GetEntries();
    Long64_t max_event = max_event_data+max_event_background;
    
    Bool_t stop_now = false;
    Long64_t max_event_testing = 12e6;
    start = std::clock();
    
    ////////////////////////Loop over Background///////////////////
    while (myReader_background.Next()) {
        
        background_loop(); //fills background histograms
        
        if (for_testing==true){
            if(timer%100000==0){
                cout << "Event #: " << timer  << " Progress: " << 100.0*timer/max_event_testing << " [%] \r";
                cout.flush();
            }
            timer++;
            if(timer==(max_event_testing/2)){
                stop_now = true;
            }
        }
        else{
            if(timer%4000000==0){
                cout << "Event #: " << timer  << " Progress: " << 100.0*timer/max_event << " [%] \r";
                cout.flush();
            }
            timer++;
            
        }
        if (stop_now==true){
            break;
        }
    }
    stop_now = false;
    
    ////////////////////// Loop over all entries of the TChain.////////////////////////////
    while (myReader.Next()) {
        first_loop();
        //gate_on_excitation_spectra=true; //set to true if no excitation gate to be used
        if(bad_particle==true||bad_gamma==true){
            //skip event
        }
        else if (gate_on_excitation_spectra==true){
            //cout<<"good event" <<endl;
            //cout<<"spam"<<endl;
            second_loop();
        }
        
        
        if (for_testing==true){
            if(timer%100000==0){
                cout << "Event #: " << timer  << " Progress: " << 100.0*timer/max_event_testing << " [%] \r";
                cout.flush();
            }
            timer++;
            if(timer==max_event_testing){
                stop_now=true;
            }
        }
        else{
            if(timer%4000000==0){
                cout << "Event #: " << timer  << " Progress: " << 100.0*timer/max_event << " [%] \r";
                cout.flush();
            }
            timer++;
        }
        if (stop_now==true){
            break;
        }
    }
    
    ////////////////write and close histogram file////////////////////
    Int_t j=0;
    Int_t i=0;
    for(j=0;j<4;j++){
        for(i=0;i<16;i++){
            hPixieADC_subtracted[j][i] -> Add(hPixieADC[j][i],1);
            //hPixieADC_subtracted[j][i] -> Add(backgroundADC[j][i],-0.11399371);//EU_before
            //hPixieADC_subtracted[j][i] -> Add(backgroundADC[j][i],-0.350194552);//EU_after
            //hPixieADC_subtracted[j][i] -> Add(backgroundADC[j][i],-2.093160377);//Ga_before
            //hPixieADC_subtracted[j][i] -> Add(backgroundADC[j][i],-2.865758754);//Ga_after
            //hPixieADC_subtracted[j][i] -> Add(backgroundADC[j][i],-0.506839622);//13C_before
            hPixieADC_subtracted[j][i] -> Add(backgroundADC[j][i],-0.667314394);//13C_after
            hPixieADC_pg[j][i]->Add(time_background[j][i],-(10.0/600.0));
        }
    }
    //cout<<gamma_amount_of_detectors_count<<endl;
    f_hist.Write();
    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    cout << "Timer: " << duration << " [s] to analyze " << max_event << " events." << endl;
    f_hist.Close();
    
}

void background_loop(){
    
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

////////////////////////first loop, checks for add-back and gates//////////////////
void first_loop(){
    
    Int_t entry_id;
    Int_t entry_size;
    Int_t sector_counter = 0;
    Int_t gamma_counter=0;
    Int_t ring_counter=0;
    Double_t time;
    Double_t en;
    Double_t calib_en;
    
    //gates
    gate_on_excitation_spectra = false;
    Double_t gamma_amount_of_detectors_check =0;
    ring_time=999999;
    sector_time=99999;
    gamma_time=99999;
    bad_gamma=false;
    
    //energy_sharing
    Double_t old_sharing_ring_en=0;
    Double_t old_sharing_ring_chn=0;
    Double_t old_sharing_ring_mod=0;
    Double_t old_sharing_sector_en=0;
    Double_t old_sharing_sector_chn=0;
    bad_particle=false;
    energy_sharing_ring=false;
    energy_sharing_sector=false;
    ring_sector_time_gate=false;
    
    //Add-back for clover 1
    Double_t ab_oldtime_cl1=0;
    Double_t ab_olden_cl1=0;
    Double_t ab_chn_cl1=0;
    Double_t count_cl1=0;
    Double_t ab_oldlocation_cl1=0;
    ab_event_cl1=false;
    
    //Add-back for clover 3
    Double_t ab_oldtime_cl3=0;
    Double_t ab_olden_cl3=0;
    Double_t ab_chn_cl3=0;
    Double_t count_cl3=0;
    Double_t ab_oldlocation_cl3=0;
    ab_event_cl3=false;
    
    for (entry_id = 0, entry_size = TRenergy.GetSize(); entry_id < entry_size; ++entry_id) {
        
        Int_t modchan=TRmodchan[entry_id];
        Int_t mod=modchan/16;
        Int_t chn=modchan%16;
        time=TRtime[entry_id];
        en=TRenergy[entry_id];
        calib_en=do_calibrate(en,mod,chn);
        //////////////////////////////////////////gates///////////////////////////////////////
        
        ///////////////energy sharing////////////////////////
        
        if( (mod==1 && (chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==2 && (chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==3 && (chn==7||chn==8||chn==9||chn==10||chn==12||chn==13||chn==14||chn==15)) ){
            if (ring_counter==0){
                ring_time=time;
                ring_mod = mod;
                ring_chn = chn;
                old_sharing_ring_en=calib_en;
                old_sharing_ring_chn=chn;
                old_sharing_ring_mod=mod;
            }
            
            if(ring_counter==1){
                //cout <<old_sharing_ring_chn-chn<<" " << old_sharing_ring_chn<< " " << chn << endl;
                if(abs(old_sharing_ring_chn-chn)==1){
                    energy_sharing_ring_en=old_sharing_sector_en+en;
                    energy_sharing_ring_chn_use=chn;
                    energy_sharing_ring_chn_ignore=old_sharing_ring_chn;
                    energy_sharing_ring_mod_use=mod;
                    energy_sharing_ring_mod_ignore=old_sharing_ring_mod;
                    energy_sharing_ring=true;
                    ring_mod = mod;
                    ring_chn = chn;
                    //cout <<old_sharing_ring_chn-chn<<" " << old_sharing_ring_chn<< " " << chn << endl;
                }
                else{
                    bad_particle=true;
                    //cout <<"bad_particle"<<endl;
                }
            }
            if(ring_counter>1){
                bad_particle=true;
                //cout <<"bad_particle"<<endl;
                
            }
            
            //////////excitation spectra gate//////////////////
            if( (calib_en>lower_lim && calib_en<upper_lim)||(energy_sharing_ring_en>lower_lim && energy_sharing_ring_en<upper_lim &&energy_sharing_ring==true)){
                gate_on_excitation_spectra=true;
            }
            if((energy_sharing_ring_en<lower_lim||energy_sharing_ring_en>upper_lim)&&(energy_sharing_ring==true)){
                bad_particle=true;
            }
            
            
            ring_counter++;
            
        }
        
        if(mod==0){
            if (sector_counter==0){
                sector_time=time;
                old_sharing_sector_en=calib_en;
                old_sharing_sector_chn=chn;
                sector_mod = mod;
                sector_chn = chn;
            }
            if (sector_counter==1){
                if(abs(old_sharing_sector_chn-chn)==1){
                    energy_sharing_sector_en=old_sharing_sector_en+en;
                    energy_sharing_sector_chn_use=chn;
                    energy_sharing_sector_chn_ignore=old_sharing_sector_chn;
                    energy_sharing_sector=true;
                    sector_mod = mod;
                    sector_chn = chn;
                }
                else{
                    bad_particle=true;
                }
            }
            if(sector_counter>1){
                bad_particle=true;
                
            }
            sector_counter++;
            
        }
        
        
        //////////////////////add-back////////////////////
        if (mod==1 && (chn==0 || chn==1 || chn==2 || chn==3) && calib_en < 10000 && calib_en>10){
            
            if ((time-ab_oldtime_cl1)<101.0 && count_cl1==1){
                ab_event_cl1=true;
                ab_use_energy_cl1=ab_olden_cl1+calib_en;
                ab_use_crystal_cl1=ab_oldlocation_cl1;
                ab_ignore_crystal_cl1=chn;
                
            }
            if(count_cl1>2||ab_olden_cl1+calib_en>10000){
                ab_event_cl1=false;
            }
            ab_oldtime_cl1=time;
            ab_olden_cl1=calib_en;
            ab_oldlocation_cl1=chn;
            count_cl1++;
        }
        
        if (mod==3 && (chn==0 || chn==1 || chn==2 || chn==3) && calib_en < 10000 && calib_en>10){
            
            if ((time-ab_oldtime_cl3)<101.0 && count_cl3==1){
                ab_event_cl3=true;
                ab_use_energy_cl3=ab_olden_cl3+calib_en;
                ab_use_crystal_cl3=ab_oldlocation_cl3;
                ab_ignore_crystal_cl3=chn;
            }
            if(count_cl3>2||ab_olden_cl3+calib_en>10000){
                ab_event_cl3=false;
            }
            ab_oldtime_cl3=time;
            ab_olden_cl3=calib_en;
            ab_oldlocation_cl3=chn;
            count_cl3++;
        }
        
        if( (mod==1 && (chn==0 || chn==1 || chn==2 || chn==3)) || (mod==3 && (chn==0 || chn==1 || chn==2 || chn==3||chn==5)) || (mod==2 && (chn==0 || chn==1 || chn==2 || chn==3 || chn==4 || chn==5)) ){
            if (gamma_counter==0){
                gamma_time=time;
            }
            if(ab_event_cl3==true||ab_event_cl1==true){
                bad_gamma=false;
            }
            else if(abs(gamma_amount_of_detectors_check-mod)>0 && gamma_counter>0){
                bad_gamma=true;
            }
            gamma_counter++;
            gamma_amount_of_detectors_check=mod;
        }
    }
    if( (ring_time-sector_time+200)>180 && (ring_time-sector_time+200)<220){
        ring_sector_time_gate=true;
    }
    
}


////////////////////////second loop fills histogram//////////////////
void second_loop(){
    
    Int_t entry_id;
    Int_t entry_size;
    Double_t time;
    Double_t en;
    Double_t calib_en;
    Double_t time_gate_lower=198;
    Double_t time_gate_upper=211;
    Double_t background_gate_lower=399;
    Double_t background_gate_upper=1001;
    Double_t fill_energy;
    Bool_t allready_filled;
    /*
     Double_t ab_use_energy_cl1;
     Double_t ab_use_crystal_cl1;
     Double_t ab_ignore_crystal_cl1;
     Double_t ab_use_energy_cl3;
     Double_t ab_use_crystal_cl3;
     Double_t ab_ignore_crystal_cl3;
     */
    
    pg_time->Fill(gamma_time-ring_time+200);
    //if (ring_sector_time_gate==true){
    ring_sector_time->Fill(sector_time-ring_time+200);
    //}
    for (entry_id = 0, entry_size = TRenergy.GetSize(); entry_id < entry_size; ++entry_id) {
        allready_filled = false;
        Int_t modchan=TRmodchan[entry_id];
        Int_t mod=modchan/16;
        Int_t chn=modchan%16;
        time=TRtime[entry_id];
        en=TRenergy[entry_id];
        calib_en=do_calibrate(en,mod,chn);
        
        if( (mod==1 && (chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==2 && (chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==3 && (chn==7||chn==8||chn==9||chn==10||chn==12||chn==13||chn==14||chn==15)) || mod==0 ){
            time_Si->Fill(time+200);
        }
        
        if( (mod==1&&(chn==0||chn==1||chn==2||chn==3))||(mod==3&&(chn==0||chn==1||chn==2||chn==3||chn==5)) ){
            // if ((time+200)>195 && (time+200)<210 && ring_sector_time_gate==true){
            time_Ge->Fill(time+200);
            //}
        }
        if( mod==2&&(chn==0||chn==1||chn==2||chn==3) ){
            //if ((time+200)>195 && (time+200)<210 && ring_sector_time_gate==true){
            time_La->Fill(time+200);
            //}
        }
        
        if( mod==2&&(chn==4||chn==5) ){
            //if ((time+200)>195 && (time+200)<210 && ring_sector_time_gate==true){
            time_La_army->Fill(time+200);
            //}
        }
        
        hPixieADC[mod][chn]->Fill(calib_en);
        //if (mod==1&&chn==13&&calib_en<1200){
        //    cout<<"spam1"<<endl;
        //}
        
        if (ab_event_cl1==true && chn==ab_use_crystal_cl1 && mod==1){
            fill_energy=do_doppler(ab_use_energy_cl1, mod, chn);
            if (abs(gamma_time-ring_time)<15 && (time+200)>time_gate_lower && (time+200)<time_gate_upper && ring_sector_time_gate==true){
                hPixieADC_pg[mod][chn]->Fill(fill_energy);
                hPixieADC_ab[mod][chn]->Fill(ab_use_energy_cl1);
            }
            if (abs(gamma_time-ring_time)<15 && (time+200)>background_gate_lower && (time+200)<background_gate_upper && ring_sector_time_gate==true){
                time_background[mod][chn]->Fill(fill_energy);
            }
            allready_filled = true;
        }
        else if(ab_event_cl1==true && chn==ab_ignore_crystal_cl1 && mod==1){
            //do nothing
            allready_filled = true;
        }
        if (ab_event_cl3==true && chn==ab_use_crystal_cl3 && mod==3 ){
            fill_energy=do_doppler(ab_use_energy_cl3, mod, chn);
            
            if (abs(gamma_time-ring_time)<15 && (time+200)>time_gate_lower && (time+200)<time_gate_upper && ring_sector_time_gate==true){
                hPixieADC_pg[mod][chn]->Fill(fill_energy);
                hPixieADC_ab[mod][chn]->Fill(ab_use_energy_cl3);
            }
            if (abs(gamma_time-ring_time)<15 && (time+200)>background_gate_lower && (time+200)<background_gate_upper && ring_sector_time_gate==true){
                time_background[mod][chn]->Fill(fill_energy);
            }
            allready_filled = true;
        }
        else if(ab_event_cl3==true && chn==ab_ignore_crystal_cl3 && mod==3){
            //do nothing
            allready_filled = true;
        }
        if (energy_sharing_sector==true && chn==energy_sharing_sector_chn_use && mod==0){
            
            if (abs(gamma_time-ring_time)<15 && (time+200)>time_gate_lower && (time+200)<time_gate_upper && ring_sector_time_gate==true){
                hPixieADC_ab[mod][chn]->Fill(energy_sharing_sector_en);
                hPixieADC_pg[mod][chn]->Fill(energy_sharing_sector_en);
            }
            if (abs(gamma_time-ring_time)<15 && (time+200)>background_gate_lower && (time+200)<background_gate_upper && ring_sector_time_gate==true){
                time_background[mod][chn]->Fill(energy_sharing_sector_en);
            }
            allready_filled = true;
        }
        else if(energy_sharing_sector==true && chn==energy_sharing_sector_chn_ignore && mod==0){
            //do nothing
            allready_filled = true;
        }
        
        if (energy_sharing_ring==true && chn==energy_sharing_ring_chn_use && mod==energy_sharing_ring_mod_use){
            
            if (abs(gamma_time-ring_time)<15 && (time+200)>time_gate_lower && (time+200)<time_gate_upper && ring_sector_time_gate==true){
                hPixieADC_pg[mod][chn]->Fill(energy_sharing_ring_en);
                hPixieADC_ab[mod][chn]->Fill(energy_sharing_ring_en);
            }
            if (abs(gamma_time-ring_time)<15 && (time+200)>background_gate_lower && (time+200)<background_gate_upper && ring_sector_time_gate==true){
                time_background[mod][chn]->Fill(energy_sharing_ring_en);
            }
            allready_filled = true;
        }
        else if(energy_sharing_ring==true && chn==energy_sharing_ring_chn_ignore && mod==energy_sharing_ring_mod_ignore){
            //do nothing
            //cout<<"spam2"<<endl;
            allready_filled = true;
        }
        if(allready_filled==false){
            fill_energy=do_doppler(calib_en, mod, chn);
            // if (calib_en>4100 && calib_en < 4600){
            //     cout<<calib_en << " " << fill_energy << endl;
            // }
            
            if (abs(gamma_time-ring_time)<15 && (time+200)>time_gate_lower && (time+200)<time_gate_upper && ring_sector_time_gate==true){
                if( (mod==1&&(chn==0||chn==1||chn==2||chn==3))||(mod==3&&(chn==0||chn==1||chn==2||chn==3||chn==5))||(mod==2&&(chn==0||chn==1||chn==2||chn==3||chn==4||chn==5)) ){
                    hPixieADC_pg[mod][chn]->Fill(fill_energy);
                    hPixieADC_ab[mod][chn]->Fill(calib_en);
                }
                else{
                    hPixieADC_pg[mod][chn]->Fill(calib_en);
                    hPixieADC_ab[mod][chn]->Fill(calib_en);
                }
            }
            if (abs(gamma_time-ring_time)<15 && (time+200)>background_gate_lower && (time+200)<background_gate_upper && ring_sector_time_gate==true){
                if( (mod==1&&(chn==0||chn==1||chn==2||chn==3))||(mod==3&&(chn==0||chn==1||chn==2||chn==3||chn==5))||(mod==2&&(chn==0||chn==1||chn==2||chn==3||chn==4||chn==5)) ){
                    time_background[mod][chn]->Fill(fill_energy);
                }
                else{
                    time_background[mod][chn]->Fill(calib_en);
                }
                
            }
            
        }
        
    }
}

void read_calibrate_file() {
    Double_t a = 0.0;
    Double_t b = 1.0;
    Double_t c = 1.0;
    Double_t d = 0.0;
    Double_t e = 0.0;
    //FILE *calib = fopen("Calibration_after_exp_13C.txt","r");
    FILE *calib = fopen("Calibration_empty.txt","r");
    for(int m=0;m<4;m++){
        for(int n=0;n<16;n++){ //The format is: offset(a) slope^1(b) slope(c)^2 phi(d) theta(e)
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
    cal_energy=(raw_energy+random_num)*(raw_energy+random_num)*slope2 + (raw_energy+random_num)*slope1 + offset;
    //cal_energy=raw_energy+random_num;
    //if (module==2&&channel==14){
    //    cout << raw_energy << " " << cal_energy << endl;
    //}
    return cal_energy;
}

/*
 Double_t do_beta(Double_t ring_en){
 Double_t beta;
 beta=sqrt(((ring_en/1000.0)*2.0)/(931.5*14.0));
 return beta;
 }
 */

Double_t do_doppler(Double_t en, Int_t gamma_mod, Int_t gamma_chn){
    Double_t dop_cor_en;
    Double_t cos_alpha;
    cos_alpha=sin(theta[ring_mod][ring_chn])*sin(theta[gamma_mod][gamma_chn])*cos(phi[sector_mod][sector_chn]-phi[gamma_mod][gamma_chn]) + cos(theta[ring_mod][ring_chn])*cos(theta[gamma_mod][gamma_chn]);
    dop_cor_en=(en)*(1.0+beta*cos_alpha)/(sqrt(1-beta*beta));
    if (gamma_mod==2&&gamma_chn==5&&en>8800&&en<9000){
        cout << "cos alpha " << cos_alpha << " energy before " << en << " doppler shifted energy " << dop_cor_en << endl;
        cout << "theta_part " << theta[ring_mod][ring_chn]  << " phi_part " << phi[sector_mod][sector_chn] << " phi_gamma " << phi[gamma_mod][gamma_chn]<< " theta_gamma " << theta[gamma_mod][gamma_chn] << endl;
        cout << " " <<endl;
        
    }
    return dop_cor_en;
}


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
    
    time_Ge = new TH1F("time_Ge","",1500,0,1500);
    time_La = new TH1F("time_La","",1500,0,1500);
    time_La_army = new TH1F("time_La_army","",1500,0,1500);
    ring_sector_time = new TH1F("time_ring_sector","",1500,0,1500);
    pg_time = new TH1F("time_particle_gamma","",1500,0,1500);
    time_Si = new TH1F("time_Si","",1500,0,1500);
    
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            hPixieADC_pg[j][i] = new TH1F(Form("ADC_complete_%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            hPixieADC_ab[j][i] = new TH1F(Form("ADC_ab%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
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
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            time_background[j][i] = new TH1F(Form("time_background%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
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
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0663-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0664-00.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0664-01.root");*/
    
    ///////////Si calibration/////////
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0819-00.root");//38 MeV
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0820-00.root");//30 MeV
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0821-00.root");//52.6MeV
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0817-00.root");//45 MeV
    
    ///////////13C after experiment/////////
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0918-00.root");
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
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/13C/run-0921-01.root");
    
    /////////////152Eu//////////////////
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0668-00.root");//58min run
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0926-00.root"); //3 hours
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0926-01.root");
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0926-02.root");
    
    ///////////background runs before//////////////
    /*Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0670-00.root");//8.48hours
     Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0670-01.root");
     Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0671-00.root");
     Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0672-00.root");
     Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0672-01.root");*/
    
    ///////////background runs after//////////////
    Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-00.root");//08:33:38
    /*Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-02.root");
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
     Background.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0917-33.root");*/
    
    
    //////////////Ga66 after experiment beamoff/////////
    /* DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0923-00.root");//10min
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0924-00.root");//1 day and 33min
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
    
    /////////////Th/////////////////
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/calibration/run-0667-00.root");
    
}
