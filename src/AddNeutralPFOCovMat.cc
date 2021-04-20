#include "AddNeutralPFOCovMat.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
#include <UTIL/LCRelationNavigator.h>
#include "EVENT/MCParticle.h"
#include "EVENT/Cluster.h"
#include "EVENT/ReconstructedParticle.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "IMPL/ParticleIDImpl.h"
#include "UTIL/PIDHandler.h"
#include "marlin/VerbosityLevels.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TTree.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;

AddFourMomentumCovMatAllPFOs aAddFourMomentumCovMatAllPFOs;

AddFourMomentumCovMatAllPFOs::AddFourMomentumCovMatAllPFOs() :

Processor("AddFourMomentumCovMatAllPFOs"),
m_nRun(0),
m_nEvt(0),
m_nRunSum(0),
m_nEvtSum(0),
m_pTFile(NULL),
m_pTTree(NULL),
m_Histograms(NULL),
m_CovMatElements(NULL),
m_NeutralPFOswithoutTrak(NULL),
m_Photon(NULL),
m_NeutralPFO(NULL),
h_clusterE_pfoE(NULL),
h_SigmaPx2(NULL),
h_SigmaPxPy(NULL),
h_SigmaPy2(NULL),
h_SigmaPxPz(NULL),
h_SigmaPyPz(NULL),
h_SigmaPz2(NULL),
h_SigmaPxE(NULL),
h_SigmaPyE(NULL),
h_SigmaPzE(NULL),
h_SigmaE2(NULL),
h_NeutPFO_PDG(NULL),
h_NeutPFO_TYPE(NULL),
h_NeutPFO_IDasPhoton(NULL),
h_NeutPFO_IDasOther(NULL),
h_NeutPFO_Weight(NULL),
h_ResidualEnergy_ph(NULL),
h_ResidualTheta_ph(NULL),
h_ResidualPhi_ph(NULL),
h_ErrorEnergy_ph(NULL),
h_ErrorTheta_ph(NULL),
h_ErrorPhi_ph(NULL),
h_NormalizedResidualEnergy_ph(NULL),
h_NormalizedResidualTheta_ph(NULL),
h_NormalizedResidualPhi_ph(NULL),
h_ResidualEnergy_NH(NULL),
h_ResidualTheta_NH(NULL),
h_ResidualPhi_NH(NULL),
h_ErrorEnergy_NH(NULL),
h_ErrorTheta_NH(NULL),
h_ErrorPhi_NH(NULL),
h_NormalizedResidualEnergy_NH(NULL),
h_NormalizedResidualTheta_NH(NULL),
h_NormalizedResidualPhi_NH(NULL),
h_NH_EclusterPlusMass_Emcp(NULL),
h_NHEnergy(NULL)
{
	_description = "Set the convariance matrix in (P,E) for all pfos (charged particles, neutral hadrons and photons)";

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"inputPfoCollection",
					"Name of input pfo collection",
					m_inputPfoCollection,
					std::string("PandoraPFOs")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"ClusterMCTruthLinkCollection",
					"Name of input m_ClusterMCTruthLink Collection",
					m_ClusterMCTruthLinkCollection,
					std::string("ClusterMCTruthLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthClusterLinkCollection",
					"Name of input MCTruthClusterLink Collection",
					m_MCTruthClusterLinkCollection,
					std::string("MCTruthClusterLink")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"outputPfoCollection",
					"Name of output pfo collection",
					m_outputPfoCollection,
					std::string("CorrectedPfoCollection")
				);

	registerProcessorParameter(	"AssumeNeutralPFOMassive",
					"true: Neutral PFOs are taken massive, false: Neutral PFOs are taken massless",
					m_AssumeNeutralPFOMassive,
					bool(true)
				);
	registerProcessorParameter(	"isClusterEnergyKinEnergy",
					"true: the cluster energy is interpreted as kinetic energy of PFO, false: the cluster energy is interpreted as momentum magnitude of PFO",
					m_isClusterEnergyKinEnergy,
					bool(false)
				);

	registerProcessorParameter(	"updatePFO4Momentum",
					"true: Update 4-momentum of PFOs, false: set 4-momentum for PFOs same as input PFO",
					m_updatePFO4Momentum,
					bool(false)
				);

	registerProcessorParameter(	"storeRootTree",
					"whether store output comparison in RootTree or not",
					m_storeRootTree,
					bool(false)
				);

	registerProcessorParameter(	"RootFile",
	                                "Name of the output root file",
					m_rootFile,
					std::string("FourMomentumCovMatAllPFOs.root")
				);

}

void AddFourMomentumCovMatAllPFOs::init()
{

	streamlog_out(MESSAGE) << "   init called  " << std::endl;
	m_Bfield = MarlinUtil::getBzAtOrigin();
	printParameters();
	streamlog_out(MESSAGE) << " BField =  "<< m_Bfield << " Tesla" << std::endl ;
	m_nRun = 0 ;
	m_nEvt = 0 ;
	m_nRunSum = 0;
	m_nEvtSum = 0;
	m_pTFile = new TFile(m_rootFile.c_str(), "recreate");
	m_pTTree = new TTree("CovMatAllPFOs", "CovMatAllPFOs");
	m_pTTree->SetDirectory(m_pTFile);
	m_pTTree->Branch("foundLinkedMCP",&m_foundLinkedMCP);
	m_pTTree->Branch("mcEnergy",&m_mcEnergy);
	m_pTTree->Branch("mcTheta",&m_mcTheta);
	m_pTTree->Branch("mcPhi",&m_mcPhi);
	m_pTTree->Branch("RecoEnergy",&m_RecoEnergy);
	m_pTTree->Branch("RecoTheta",&m_RecoTheta);
	m_pTTree->Branch("RecoPhi",&m_RecoPhi);
	m_pTTree->Branch("ResidualEnergy",&m_ResidualEnergy);
	m_pTTree->Branch("ResidualTheta",&m_ResidualTheta);
	m_pTTree->Branch("ResidualPhi",&m_ResidualPhi);
	m_pTTree->Branch("ErrorEnergy",&m_ErrorEnergy);
	m_pTTree->Branch("ErrorTheta",&m_ErrorTheta);
	m_pTTree->Branch("ErrorPhi",&m_ErrorPhi);
	m_pTTree->Branch("NormalizedResidualEnergy",&m_NormalizedResidualEnergy);
	m_pTTree->Branch("NormalizedResidualTheta",&m_NormalizedResidualTheta);
	m_pTTree->Branch("NormalizedResidualPhi",&m_NormalizedResidualPhi);
	m_pTTree->Branch("foundLinkedMCP_Ph",&m_foundLinkedMCP_Ph);
	m_pTTree->Branch("mcEnergy_Ph",&m_mcEnergy_Ph);
	m_pTTree->Branch("mcTheta_Ph",&m_mcTheta_Ph);
	m_pTTree->Branch("mcPhi_Ph",&m_mcPhi_Ph);
	m_pTTree->Branch("RecoEnergy_Ph",&m_RecoEnergy_Ph);
	m_pTTree->Branch("RecoTheta_Ph",&m_RecoTheta_Ph);
	m_pTTree->Branch("RecoPhi_Ph",&m_RecoPhi_Ph);
	m_pTTree->Branch("ResidualEnergy_Ph",&m_ResidualEnergy_Ph);
	m_pTTree->Branch("ResidualTheta_Ph",&m_ResidualTheta_Ph);
	m_pTTree->Branch("ResidualPhi_Ph",&m_ResidualPhi_Ph);
	m_pTTree->Branch("ErrorEnergy_Ph",&m_ErrorEnergy_Ph);
	m_pTTree->Branch("ErrorTheta_Ph",&m_ErrorTheta_Ph);
	m_pTTree->Branch("ErrorPhi_Ph",&m_ErrorPhi_Ph);
	m_pTTree->Branch("NormalizedResidualEnergy_Ph",&m_NormalizedResidualEnergy_Ph);
	m_pTTree->Branch("NormalizedResidualTheta_Ph",&m_NormalizedResidualTheta_Ph);
	m_pTTree->Branch("NormalizedResidualPhi_Ph",&m_NormalizedResidualPhi_Ph);
	m_pTTree->Branch("foundLinkedMCP_NH",&m_foundLinkedMCP_NH);
	m_pTTree->Branch("mcEnergy_NH",&m_mcEnergy_NH);
	m_pTTree->Branch("mcTheta_NH",&m_mcTheta_NH);
	m_pTTree->Branch("mcPhi_NH",&m_mcPhi_NH);
	m_pTTree->Branch("RecoEnergy_NH",&m_RecoEnergy_NH);
	m_pTTree->Branch("RecoTheta_NH",&m_RecoTheta_NH);
	m_pTTree->Branch("RecoPhi_NH",&m_RecoPhi_NH);
	m_pTTree->Branch("ResidualEnergy_NH",&m_ResidualEnergy_NH);
	m_pTTree->Branch("ResidualTheta_NH",&m_ResidualTheta_NH);
	m_pTTree->Branch("ResidualPhi_NH",&m_ResidualPhi_NH);
	m_pTTree->Branch("ErrorEnergy_NH",&m_ErrorEnergy_NH);
	m_pTTree->Branch("ErrorTheta_NH",&m_ErrorTheta_NH);
	m_pTTree->Branch("ErrorPhi_NH",&m_ErrorPhi_NH);
	m_pTTree->Branch("NormalizedResidualEnergy_NH",&m_NormalizedResidualEnergy_NH);
	m_pTTree->Branch("NormalizedResidualTheta_NH",&m_NormalizedResidualTheta_NH);
	m_pTTree->Branch("NormalizedResidualPhi_NH",&m_NormalizedResidualPhi_NH);
	h_clusterE_pfoE = new TH2F("Neutral PFOs (n_{Tracks} = 0)", "; E_{Cluster} [GeV]; E_{PFO} [GeV]", 10000, 0.0, 100., 10000, 0.0, 100.);
	h_SigmaPx2 = new TH2F("Charged PFOs", "; #sigma_{p_{x}}^{2} (new PFO) [GeV^{2}]; #sigma_{p_{x}}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPxPy = new TH2F("Charged PFOs", "; #sigma_{p_{x}p_{y}} (new PFO) [GeV^{2}]; #sigma_{p_{x}p_{y}} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPy2 = new TH2F("Charged PFOs", "; #sigma_{p_{y}}^{2} (new PFO) [GeV^{2}]; #sigma_{p_{y}}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPxPz = new TH2F("Charged PFOs", "; #sigma_{p_{x}p_{z}} (new PFO) [GeV^{2}]; #sigma_{p_{x}p_{z}} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyPz = new TH2F("Charged PFOs", "; #sigma_{p_{y}p_{z}} (new PFO) [GeV^{2}]; #sigma_{p_{y}p_{z}} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPz2 = new TH2F("Charged PFOs", "; #sigma_{p_{z}}^{2} (new PFO) [GeV^{2}]; #sigma_{p_{z}}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.1, 400, 0.0, 0.1);
	h_SigmaPxE = new TH2F("Charged PFOs", "; #sigma_{p_{x}E} (new PFO) [GeV^{2}]; #sigma_{p_{x}E} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPyE = new TH2F("Charged PFOs", "; #sigma_{p_{y}E} (new PFO) [GeV^{2}]; #sigma_{p_{y}E} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaPzE = new TH2F("Charged PFOs", "; #sigma_{p_{z}E} (new PFO) [GeV^{2}]; #sigma_{p_{z}E} (old PFO) [GeV^{2}]", 400, -0.1, 0.1, 400, -0.1, 0.1);
	h_SigmaE2 = new TH2F("Charged PFOs", "; #sigma_{E}^{2} (new PFO) [GeV^{2}]; #sigma_{E}^{2} (old PFO) [GeV^{2}]", 400, 0.0, 0.6, 400, 0.0, 0.6);
	h_NeutPFO_PDG = new TH1I("Neutral PFOs PDG", "; PDG Code", 200001, -100000.5, 100000.5);
	h_NeutPFO_TYPE = new TH1I("Neutral PFOs TYPE", "; True Part. Type", 15, 0, 15);
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(1,"e^{#pm}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(3,"#gamma");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(4,"K^{0}_{L}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(5,"#pi^{#pm}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(6,"K^{0}_{S}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(7,"K^{#pm}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(8,"n");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(9,"p");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(10,"#Sigma^{-}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(11,"#Lambda");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(12,"#Sigma^{+}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(13,"#Xi^{-}");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(14,"#Xi");
	h_NeutPFO_TYPE->GetXaxis()->SetBinLabel(15,"Others");
	h_NeutPFO_IDasPhoton = new TH1I("Photons", "; True Part. Type", 15, 0, 15);
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(1,"e^{#pm}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(3,"#gamma");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(4,"K^{0}_{L}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(5,"#pi^{#pm}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(6,"K^{0}_{S}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(7,"K^{#pm}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(8,"n");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(9,"p");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(10,"#Sigma^{-}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(11,"#Lambda");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(12,"#Sigma^{+}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(13,"#Xi^{-}");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(14,"#Xi");
	h_NeutPFO_IDasPhoton->GetXaxis()->SetBinLabel(15,"Others");
	h_NeutPFO_IDasOther = new TH1I("Other Neutal PFOs", "; True Part. Type", 15, 0, 15);
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(1,"e^{#pm}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(3,"#gamma");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(4,"K^{0}_{L}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(5,"#pi^{#pm}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(6,"K^{0}_{S}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(7,"K^{#pm}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(8,"n");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(9,"p");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(10,"#Sigma^{-}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(11,"#Lambda");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(12,"#Sigma^{+}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(13,"#Xi^{-}");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(14,"#Xi");
	h_NeutPFO_IDasOther->GetXaxis()->SetBinLabel(15,"Others");
	h_NeutPFO_Weight = new TH1F("Neutral Hadrons MCP Link Weight", "; Link weight", 100, 0.0, 1.0);
	h_ResidualEnergy_ph = new TH1F("Photons", "; E_{REC} - E_{MCP} [GeV]", 200, -10.0, 10.0);
	h_ResidualTheta_ph = new TH1F("Photons", "; #theta_{REC} - #theta_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ResidualPhi_ph = new TH1F("Photons", "; #phi_{REC} - #phi_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ErrorEnergy_ph = new TH1F("Photons", "; #sigma_{E} [GeV]", 1000, 0.0, 10.0);
	h_ErrorTheta_ph = new TH1F("Photons", "; #sigma_{#theta} [radian]", 10000, 0.0, 1.0);
	h_ErrorPhi_ph = new TH1F("Photons", "; #sigma_{#phi} [radian]", 10000, 0.0, 1.0);
	h_NormalizedResidualEnergy_ph = new TH1F("Photons", "; (E_{REC} - E_{MCP}) / #sigma_{E}", 200, -10.0, 10.0);
	h_NormalizedResidualTheta_ph = new TH1F("Photons", "; (#theta_{REC} - #theta_{MCP}) / #sigma_{#theta}", 200, -10.0, 10.0);
	h_NormalizedResidualPhi_ph = new TH1F("Photons", "; (#phi_{REC} - #phi_{MCP}) / #sigma_{#phi}", 200, -10.0, 10.0);
	h_ResidualEnergy_NH = new TH1F("Neutral Hadrons", "; E_{REC} - E_{MCP} [GeV]", 200, -10.0, 10.0);
	h_ResidualTheta_NH = new TH1F("Neutral Hadrons", "; #theta_{REC} - #theta_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ResidualPhi_NH = new TH1F("Neutral Hadrons", "; #phi_{REC} - #phi_{MCP} [radian]", 6800, -3.4, 3.4);
	h_ErrorEnergy_NH = new TH1F("Neutral Hadrons", "; #sigma_{E} [GeV]", 1000, 0.0, 10.0);
	h_ErrorTheta_NH = new TH1F("Neutral Hadrons", "; #sigma_{#theta} [radian]", 10000, 0.0, 1.0);
	h_ErrorPhi_NH = new TH1F("Neutral Hadrons", "; #sigma_{#phi} [radian]", 10000, 0.0, 1.0);
	h_NormalizedResidualEnergy_NH = new TH1F("Neutral Hadrons", "; (E_{REC} - E_{MCP}) / #sigma_{E}", 200, -10.0, 10.0);
	h_NormalizedResidualTheta_NH = new TH1F("Neutral Hadrons", "; (#theta_{REC} - #theta_{MCP}) / #sigma_{#theta}", 200, -10.0, 10.0);
	h_NormalizedResidualPhi_NH = new TH1F("Neutral Hadrons", "; (#phi_{REC} - #phi_{MCP}) / #sigma_{#phi}", 200, -10.0, 10.0);
	h_NH_EclusterPlusMass_Emcp = new TH2F("Neutral Hadrons", "; E_{MCP} [GeV]; E_{cluster} + m [GeV]", 1000, 0.0, 10.0, 1000, 0.0, 10.0);
	h_NHEnergy = new TH2F("Neutral Hadrons", "; E_{MCP} [GeV]; E_{PFO} [GeV]", 1000, 0.0, 10.0, 1000, 0.0, 10.0);

}

void AddFourMomentumCovMatAllPFOs::Clear()
{
	m_foundLinkedMCP.clear();
	m_mcEnergy.clear();
	m_mcTheta.clear();
	m_mcPhi.clear();
	m_RecoEnergy.clear();
	m_RecoTheta.clear();
	m_RecoPhi.clear();
	m_ResidualEnergy.clear();
	m_ResidualTheta.clear();
	m_ResidualPhi.clear();
	m_ErrorEnergy.clear();
	m_ErrorTheta.clear();
	m_ErrorPhi.clear();
	m_NormalizedResidualEnergy.clear();
	m_NormalizedResidualTheta.clear();
	m_NormalizedResidualPhi.clear();
	m_foundLinkedMCP_Ph.clear();
	m_mcEnergy_Ph.clear();
	m_mcTheta_Ph.clear();
	m_mcPhi_Ph.clear();
	m_RecoEnergy_Ph.clear();
	m_RecoTheta_Ph.clear();
	m_RecoPhi_Ph.clear();
	m_ResidualEnergy_Ph.clear();
	m_ResidualTheta_Ph.clear();
	m_ResidualPhi_Ph.clear();
	m_ErrorEnergy_Ph.clear();
	m_ErrorTheta_Ph.clear();
	m_ErrorPhi_Ph.clear();
	m_NormalizedResidualEnergy_Ph.clear();
	m_NormalizedResidualTheta_Ph.clear();
	m_NormalizedResidualPhi_Ph.clear();
	m_foundLinkedMCP_NH.clear();
	m_mcEnergy_NH.clear();
	m_mcTheta_NH.clear();
	m_mcPhi_NH.clear();
	m_RecoEnergy_NH.clear();
	m_RecoTheta_NH.clear();
	m_RecoPhi_NH.clear();
	m_ResidualEnergy_NH.clear();
	m_ResidualTheta_NH.clear();
	m_ResidualPhi_NH.clear();
	m_ErrorEnergy_NH.clear();
	m_ErrorTheta_NH.clear();
	m_ErrorPhi_NH.clear();
	m_NormalizedResidualEnergy_NH.clear();
	m_NormalizedResidualTheta_NH.clear();
	m_NormalizedResidualPhi_NH.clear();
}

void AddFourMomentumCovMatAllPFOs::processRunHeader()
{

	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;

}

void AddFourMomentumCovMatAllPFOs::processEvent( EVENT::LCEvent *pLCEvent )
{

	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	++m_nEvtSum;
	streamlog_out(MESSAGE) << "processed event 	" << m_nEvtSum << std::endl;

	LCCollection *inputPfoCollection{};
	this->Clear();

	try
	{
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		int n_PFO = inputPfoCollection->getNumberOfElements();
		LCCollectionVec *m_col_outputPfo = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

		streamlog_out(DEBUG) << "Investigated All PFOs" << std::endl;
		m_pTTree->Fill();
	}
        catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input collection not found in event " << m_nEvt << std::endl;
        }

}

std::vector<float> AddFourMomentumCovMatAllPFOs::UpdateNeutralPFOCovMat( TVector3 clusterPosition , float pfoEc , float pfoMass , std::vector<float> clusterPositionError , float clusterEnergyError )
{

//	Obtain covariance matrix on (px,py,pz,E) from the
//	covariance matrix on cluster parameters (px,py,pz,Ek=Ec=E-E0).
//	=> E = Ek + m	;	|p| = sqrt( Ek^2 + 2mEk )
//	define the jacobian as the 4x4 matrix:
//
//
//
//			Dpx/Dx			Dpy/Dx			Dpz/Dx			DE/Dx
//
//			Dpx/Dy			Dpy/Dy			Dpz/Dy			DE/Dy
//	J =
//			Dpx/Dz			Dpy/Dz			Dpz/Dz			DE/Dz
//
//			Dpx/DEk			Dpy/DEk			Dpz/DEk			DE/DEk
//
//
//
//
//
//			 |P|.(r2-x2)/r3		-|P|.x.y/r3		-|P|.x.z/r3		0
//
//			-|P|.y.x/r3		 |P|.(r2-y2)/r3		-|P|.y.z/r3		0
//	J =
//			-|P|.z.x/r3		-|P|.z.y/r3		 |P|.(r2-z2)/r3		0
//
//			 (E/|p|).(x/r)		 (E/|p|).(y/r)		 (E/|p|).(z/r)		1
//
//
//
//
//	CovMatrix elements in terms of cluster position error and cluster energy error:
//
//			x.x			x.y			x.z			x.Ec
//
//			y.x			y.y			y.z			y.Ec
//	Cov =
//			z.x			z.y			z.z			z.Ec
//
//			Ec.x			Ec.y			Ec.z			Ec.Ec
//
//
//

	const int rows			= 4; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int kspace_time_dim	= 4;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	std::vector<float> covP;

//	pfoMass			= 0.0;

	float pfoX		=	clusterPosition.X();
	float pfoY		=	clusterPosition.Y();
	float pfoZ		=	clusterPosition.Z();
	float pfoR		=	std::sqrt( pow( pfoX , 2 ) + pow( pfoY , 2 ) + pow( pfoZ , 2 ) );
	float pfoX2		=	pow( pfoX , 2 );
	float pfoY2		=	pow( pfoY , 2 );
	float pfoZ2		=	pow( pfoZ , 2 );
	float pfoR2		=	pow( pfoR , 2 );
	float pfoR3		=	pow( pfoR , 3 );
	float SigmaX2		=	clusterPositionError[ 0 ];
	float SigmaXY		=	clusterPositionError[ 1 ];
	float SigmaY2		=	clusterPositionError[ 2 ];
	float SigmaXZ		=	clusterPositionError[ 3 ];
	float SigmaYZ		=	clusterPositionError[ 4 ];
	float SigmaZ2		=	clusterPositionError[ 5 ];
	float SigmaE2		=	pow( clusterEnergyError , 2 );

	float pfoP;
	float pfoE;
	float derivative_coeff	= 1.0;
	if ( m_isClusterEnergyKinEnergy )
	{
		pfoP			= sqrt( pow( pfoEc , 2 ) + 2 * pfoMass * pfoEc );
		pfoE			= sqrt( pow( pfoP , 2 ) + pow( pfoMass , 2 ) );
//		derivative_coeff	= 1.;
	}
	else
	{
		pfoP			= pfoEc;
		pfoE			= sqrt( pow( pfoP , 2 ) + pow( pfoMass , 2 ) );
		if ( m_useTrueJacobian ) derivative_coeff	= pfoP / pfoE;
	}

	streamlog_out(DEBUG) << "Cluster information obtained" << std::endl;

//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
	{
		pfoP * ( pfoR2 - pfoX2 ) / pfoR3			,	-pfoP * pfoX * pfoY / pfoR3				,	-pfoP * pfoX * pfoZ / pfoR3				,	0			,
		-pfoP * pfoY * pfoX / pfoR3				,	pfoP * ( pfoR2 - pfoY2 ) / pfoR3			,	-pfoP * pfoY * pfoZ / pfoR3				,	0			,
		-pfoP * pfoZ * pfoX / pfoR3				,	-pfoP * pfoZ * pfoY / pfoR3				,	pfoP * ( pfoR2 - pfoZ2 ) / pfoR3			,	0			,
		derivative_coeff * pfoE * pfoX / ( pfoP * pfoR )	,	derivative_coeff * pfoE * pfoY / ( pfoP * pfoR )	,	derivative_coeff * pfoE * pfoZ / ( pfoP * pfoR )	,	derivative_coeff
	};
/*
	streamlog_out(MESSAGE) << "******************************************************************************************" << std::endl;
	streamlog_out(DEBUG) << "Jacobian array formed by rows" << std::endl;
	streamlog_out(MESSAGE) << "Jacobain CovMat(x,y,z,E) -> CovMat (Px,Py,Pz,E):" << std::endl;
	streamlog_out(MESSAGE) << "{" << std::endl;
	streamlog_out(MESSAGE) << "	" << jacobian_by_rows[ 0 ] << "	,	" << jacobian_by_rows[ 1 ] << "	,	" << jacobian_by_rows[ 2 ] << "	,	" << jacobian_by_rows[ 3 ] << std::endl;
	streamlog_out(MESSAGE) << "	" << jacobian_by_rows[ 4 ] << "	,	" << jacobian_by_rows[ 5 ] << "	,	" << jacobian_by_rows[ 6 ] << "	,	" << jacobian_by_rows[ 7 ] << std::endl;
	streamlog_out(MESSAGE) << "	" << jacobian_by_rows[ 8 ] << "	,	" << jacobian_by_rows[ 9 ] << "	,	" << jacobian_by_rows[ 10 ] << "	,	" << jacobian_by_rows[ 11 ] << std::endl;
	streamlog_out(MESSAGE) << "	" << jacobian_by_rows[ 12 ] << "	,	" << jacobian_by_rows[ 13 ] << "	,	" << jacobian_by_rows[ 14 ] << "	,	" << jacobian_by_rows[ 15 ] << std::endl;
	streamlog_out(MESSAGE) << "}" << std::endl;
	streamlog_out(MESSAGE) << "******************************************************************************************" << std::endl;
*/

//	construct the Jacobian using previous array ("F" if filling by columns, "C" if filling by rows, $ROOTSYS/math/matrix/src/TMatrixT.cxx)
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG) << "Jacobian array converted to Jacobian matrix" << std::endl;

//	cluster covariance matrix by rows
	double cluster_cov_matrix_by_rows[rows*rows] =
			{
				SigmaX2		,	SigmaXY		,	SigmaXZ		,	0	,
				SigmaXY		,	SigmaY2		,	SigmaYZ		,	0	,
				SigmaXZ		,	SigmaYZ		,	SigmaZ2		,	0	,
				0		,	0		,	0		,	SigmaE2
			};
/*
	streamlog_out(DEBUG) << "cluster covariance matrix array formed by rows" << std::endl;
	streamlog_out(MESSAGE) << "CovMat(x,y,z,E):" << std::endl;
	streamlog_out(MESSAGE) << "{" << std::endl;
	streamlog_out(MESSAGE) << "	" << cluster_cov_matrix_by_rows[ 0 ] << "	,	" << cluster_cov_matrix_by_rows[ 1 ] << "	,	" << cluster_cov_matrix_by_rows[ 2 ] << "	,	" << cluster_cov_matrix_by_rows[ 3 ] << std::endl;
	streamlog_out(MESSAGE) << "	" << cluster_cov_matrix_by_rows[ 4 ] << "	,	" << cluster_cov_matrix_by_rows[ 5 ] << "	,	" << cluster_cov_matrix_by_rows[ 6 ] << "	,	" << cluster_cov_matrix_by_rows[ 7 ] << std::endl;
	streamlog_out(MESSAGE) << "	" << cluster_cov_matrix_by_rows[ 8 ] << "	,	" << cluster_cov_matrix_by_rows[ 9 ] << "	,	" << cluster_cov_matrix_by_rows[ 10 ] << "	,	" << cluster_cov_matrix_by_rows[ 11 ] << std::endl;
	streamlog_out(MESSAGE) << "	" << cluster_cov_matrix_by_rows[ 12 ] << "	,	" << cluster_cov_matrix_by_rows[ 13 ] << "	,	" << cluster_cov_matrix_by_rows[ 14 ] << "	,	" << cluster_cov_matrix_by_rows[ 15 ] << std::endl;
	streamlog_out(MESSAGE) << "}" << std::endl;
	streamlog_out(MESSAGE) << "******************************************************************************************" << std::endl;
*/
	TMatrixD covMatrix_cluster(rows,rows, cluster_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG) << "cluster covariance matrix array converted to cluster covariance matrix" << std::endl;

	covMatrixMomenta.Mult( TMatrixD( jacobian ,
					TMatrixD::kTransposeMult ,
					covMatrix_cluster) ,
					jacobian
					);
/*
	streamlog_out(DEBUG) << "cluster covariance matrix array converted to FourMomentumCovariance matrix" << std::endl;
	streamlog_out(MESSAGE) << "CovMat(x,y,z,E):" << std::endl;
	streamlog_out(MESSAGE) << "{" << std::endl;
	streamlog_out(MESSAGE) << "	" << covMatrixMomenta(0,0) << "	,	" << covMatrixMomenta(0,1) << "	,	" << covMatrixMomenta(0,2) << "	,	" << covMatrixMomenta(0,3) << std::endl;
	streamlog_out(MESSAGE) << "	" << covMatrixMomenta(1,0) << "	,	" << covMatrixMomenta(1,1) << "	,	" << covMatrixMomenta(1,2) << "	,	" << covMatrixMomenta(1,3) << std::endl;
	streamlog_out(MESSAGE) << "	" << covMatrixMomenta(2,0) << "	,	" << covMatrixMomenta(2,1) << "	,	" << covMatrixMomenta(2,2) << "	,	" << covMatrixMomenta(2,3) << std::endl;
	streamlog_out(MESSAGE) << "	" << covMatrixMomenta(3,0) << "	,	" << covMatrixMomenta(3,1) << "	,	" << covMatrixMomenta(3,2) << "	,	" << covMatrixMomenta(3,3) << std::endl;
	streamlog_out(MESSAGE) << "}" << std::endl;
	streamlog_out(MESSAGE) << "******************************************************************************************" << std::endl;
*/
	covP.push_back( covMatrixMomenta(0,0) ); // x-x
	covP.push_back( covMatrixMomenta(1,0) ); // y-x
	covP.push_back( covMatrixMomenta(1,1) ); // y-y
	covP.push_back( covMatrixMomenta(2,0) ); // z-x
	covP.push_back( covMatrixMomenta(2,1) ); // z-y
	covP.push_back( covMatrixMomenta(2,2) ); // z-z
	covP.push_back( covMatrixMomenta(3,0) ); // e-x
	covP.push_back( covMatrixMomenta(3,1) ); // e-y
	covP.push_back( covMatrixMomenta(3,2) ); // e-z
	covP.push_back( covMatrixMomenta(3,3) ); // e-e
	streamlog_out(DEBUG) << "FourMomentumCovarianceMatrix Filled succesfully" << std::endl;

	return covP;

}

std::vector<float> AddFourMomentumCovMatAllPFOs::getClusterDirectionError( TVector3 clusterPosition , std::vector<float> clusterPositionError )
{

//	Obtain covariance matrix on (R,Theta,Phi) from the
//	covariance matrix on cluster parameters (x,y,z).
//	=> R^2 = x^2 + y^2 + z^2	;	tan(Theta) = sqrt( x^2 + y^2 ) / z	;	tan(Phi) = y / x
//	define the jacobian as the 4x4 matrix:
//
//
//
//			DR/Dx			DTheta/Dx		DPhi/Dx
//
//	J =		DR/Dy			DTheta/Dy		DPhi/Dy
//
//			DR/Dz			DTheta/Dz		DPhi/Dz
//
//
//
//
//			x/R			x.z/(R3.sqrt(1-(z2/R2)))		-y/(x2+y2)
//
//	J =		y/R			y.z/(R3.sqrt(1-(z2/R2)))		x/(x2+y2)
//
//			z/R			-sqrt(1-(z2/R2))/R			0
//
//
//
//
//
//	CovMatrix elements in terms of cluster position error:
//
//			x.x			x.y			x.z
//
//	Cov =		y.x			y.y			y.z
//
//			z.x			z.y			z.z
//
//
//
//

	const int rows			= 3; // n rows jacobian
	const int columns		= 3; // n columns jacobian
	const int kspace_time_dim	= 3;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	std::vector<float> covR;

//	pfoMass			= 0.0;

	float pfoX		=	clusterPosition.X();
	float pfoY		=	clusterPosition.Y();
	float pfoZ		=	clusterPosition.Z();
	float pfoR		=	std::sqrt( pow( pfoX , 2 ) + pow( pfoY , 2 ) + pow( pfoZ , 2 ) );
	float pfoR2		=	pow( pfoR , 2 );
	float pfoR3		=	pow( pfoR , 3 );
	float pfoX2		=	pow( pfoX , 2 );
	float pfoY2		=	pow( pfoY , 2 );
	float pfoZ2		=	pow( pfoZ , 2 );
	float SigmaX2		=	clusterPositionError[ 0 ];
	float SigmaXY		=	clusterPositionError[ 1 ];
	float SigmaY2		=	clusterPositionError[ 2 ];
	float SigmaXZ		=	clusterPositionError[ 3 ];
	float SigmaYZ		=	clusterPositionError[ 4 ];
	float SigmaZ2		=	clusterPositionError[ 5 ];

	streamlog_out(DEBUG) << "Cluster information obtained" << std::endl;

//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
			{
				pfoX / pfoR		,		pfoX * pfoZ / ( pfoR3 * sqrt( 1 - pfoZ2 / pfoR2 ) )		,	-1 * pfoY / ( pfoY2 + pfoX2 )		,
				pfoY / pfoR		,		pfoY * pfoZ / ( pfoR3 * sqrt( 1 - pfoZ2 / pfoR2 ) )		,	pfoX / ( pfoY2 + pfoX2 )		,
				pfoZ / pfoR		,		-1 * sqrt( 1 - pfoZ2 / pfoR2 ) / pfoR				,		0
			};

	streamlog_out(DEBUG) << "Jacobian array formed by rows" << std::endl;

//	construct the Jacobian using previous array ("F" if filling by columns, "C", if filling by rows, $ROOTSYS/math/matrix/src/TMatrixT.cxx)
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG) << "Jacobian array converted to Jacobian matrix" << std::endl;

//	cluster covariance matrix by rows
	double cluster_cov_matrix_by_rows[rows*rows] =
			{
				SigmaX2		,	SigmaXY		,	SigmaXZ		,
				SigmaXY		,	SigmaY2		,	SigmaYZ		,
				SigmaXZ		,	SigmaYZ		,	SigmaZ2
			};
	streamlog_out(DEBUG) << "cluster covariance matrix array formed by rows" << std::endl;

	TMatrixD covMatrix_cluster(rows,rows, cluster_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG) << "cluster covariance matrix array converted to cluster covariance matrix" << std::endl;

	covMatrixMomenta.Mult( TMatrixD( jacobian ,
					TMatrixD::kTransposeMult ,
					covMatrix_cluster) ,
					jacobian
					);
	streamlog_out(DEBUG) << "cluster covariance matrix array in cartesian coordinate system converted to cluster covariance matrix array in spherical (polar) coordinate system" << std::endl;

	covR.push_back( covMatrixMomenta(0,0) ); // R-R
	covR.push_back( covMatrixMomenta(1,0) ); // Theta-R
	covR.push_back( covMatrixMomenta(1,1) ); // Theta-Theta
	covR.push_back( covMatrixMomenta(2,0) ); // Phi-R
	covR.push_back( covMatrixMomenta(2,1) ); // Phi-Theta
	covR.push_back( covMatrixMomenta(2,2) ); // Phi-Phi
	streamlog_out(DEBUG) << "FourMomentumCovarianceMatrix Filled succesfully" << std::endl;

	return covR;

}

std::vector<float> AddFourMomentumCovMatAllPFOs::getPFOResidual( TLorentzVector pfoFourMomentum , TLorentzVector mcpFourMomentum )
{
	std::vector<float> pfoResidual;

	float pfoPx		= pfoFourMomentum.Px();
	float pfoPy		= pfoFourMomentum.Py();
	float pfoPz		= pfoFourMomentum.Pz();
	float pfoE		= pfoFourMomentum.E();
	TVector3 pfoPvec( pfoPx , pfoPy , pfoPz ); pfoPvec.SetMag(1.0);
	TVector3 pfoPTvec( pfoPx , pfoPy , 0.0 ); pfoPTvec.SetMag(1.0);
	float pfoTheta		= pfoPvec.Theta();
	float pfoPhi		= pfoPvec.Phi();
	streamlog_out(DEBUG) << "PFO 4-Momentum ( px , py , pz , E , Theta , Phi ) = ( " << pfoPx << " , " << pfoPy << " , " << pfoPz << " , " << pfoE << " , " << pfoTheta << " , " << pfoPhi << " )" << std::endl;

	float mcpPx		= mcpFourMomentum.Px();
	float mcpPy		= mcpFourMomentum.Py();
	float mcpPz		= mcpFourMomentum.Pz();
	float mcpE		= mcpFourMomentum.E();
	TVector3 mcpPvec( mcpPx , mcpPy , mcpPz ); mcpPvec.SetMag(1.0);
	TVector3 mcpPTvec( mcpPx , mcpPy , 0.0 ); mcpPTvec.SetMag(1.0);
	float mcpTheta		= mcpPvec.Theta();
	float mcpPhi		= mcpPvec.Phi();
	streamlog_out(DEBUG) << "MCP 4-Momentum ( px , py , pz , E , Theta , Phi ) = ( " << mcpPx << " , " << mcpPy << " , " << mcpPz << " , " << mcpE << " , " << mcpTheta << " , " << mcpPhi << " )" << std::endl;

	float ResidualEnergy	= pfoE - mcpE;
	float ResidualTheta	= pfoTheta - mcpTheta;
	float ResidualPhi	= 0.0;
	if ( pfoPhi > mcpPhi )
	{
		ResidualPhi	= acos( pfoPTvec.Dot( mcpPTvec ) );
	}
	else
	{
		ResidualPhi	= -acos( pfoPTvec.Dot( mcpPTvec ) );
	}
	streamlog_out(DEBUG) << "	Residuals	( deltaE , deltaTheta , deltaPhi ) = ( " << ResidualEnergy << "	,	" << ResidualTheta << "	, " << ResidualPhi << "	)" << std::endl;

	pfoResidual.push_back( ResidualEnergy );
	pfoResidual.push_back( ResidualTheta );
	pfoResidual.push_back( ResidualPhi );

	return pfoResidual;

}

std::vector<float> AddFourMomentumCovMatAllPFOs::getPFOCovMatPolarCoordinate( TLorentzVector pfoFourMomentum , std::vector<float> pfoCovMat )
{

//	Obtain covariance matrix on (Theta,Phi,P,E) from the
//	covariance matrix on (Px,Py,Pz,E).
//	=> P^2 = Px^2 + Py^2 + Pz^2	;	tan(Theta) = sqrt( Px^2 + Py^2 ) / Pz	;	tan(Phi) = Py / Px
//	define the jacobian as the 4x4 matrix:
//
//
//
//			DTheta/DPx			DPhi/DPx		DP/DPx			DE/DPx
//
//			DTheta/DPy			DPhi/DPy		DP/DPy			DE/DPy
//	J =
//			DTheta/DPz			DPhi/DPz		DP/DPz			DE/DPz
//
//			DTheta/DE			DPhi/DE		DP/DE			DE/DE
//
//
//
//
//			(Px.Pz)/(Pt.P2)		-Py/Pt2		Px/P			0
//
//			(Py.Pz)/(Pt.P2)		Px/Pt2			Py/P			0
//	J =
//			-Pt/P2				0			Pz/P			0
//
//			0				0			0			1
//
//
//
//
//
//	CovMatrix elements in terms of 4-momentum:
//
//			Px.Px			Px.Py			Px.Pz			Px.E
//
//			Py.Px			Py.Py			Py.Pz			Py.E
//	Cov =
//			Pz.Px			Pz.Py			Pz.Pz			Pz.E
//
//			E.Px			E.Py			E.Pz			E.E
//
//
//
//

	const int rows			= 4; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int kspace_time_dim	= 4;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	std::vector<float> covP;

//	pfoMass			= 0.0;

	float totPx		=	pfoFourMomentum.Px();
	float totPy		=	pfoFourMomentum.Py();
	float totPz		=	pfoFourMomentum.Pz();
	float totP		=	std::sqrt( pow( totPx , 2 ) + pow( totPy , 2 ) + pow( totPz , 2 ) );
	float totP2		=	pow( totP , 2 );
	float totPt		=	std::sqrt( pow( totPx , 2 ) + pow( totPy , 2 ) );
	float totPt2		=	pow( totPt , 2 );
	float SigmaPx2		=	pfoCovMat[ 0 ];
	float SigmaPxPy	=	pfoCovMat[ 1 ];
	float SigmaPy2		=	pfoCovMat[ 2 ];
	float SigmaPxPz	=	pfoCovMat[ 3 ];
	float SigmaPyPz	=	pfoCovMat[ 4 ];
	float SigmaPz2		=	pfoCovMat[ 5 ];
	float SigmaPxE		=	pfoCovMat[ 6 ];
	float SigmaPyE		=	pfoCovMat[ 7 ];
	float SigmaPzE		=	pfoCovMat[ 8 ];
	float SigmaE2		=	pfoCovMat[ 9 ];

	streamlog_out(DEBUG) << "Cluster information obtained" << std::endl;

//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
	{
		( totPx * totPz ) / ( totPt * totP2 )		,	-totPy / totPt2		,	totPx / totP		,	0	,
		( totPy * totPz ) / ( totPt * totP2 )		,	totPx / totPt2		,	totPy / totP		,	0	,
		-totPt / totP2					,	0			,	totPz / totP		,	0	,
		0						,	0			,	0			,	1
	};

	streamlog_out(DEBUG) << "Jacobian array formed by rows" << std::endl;

//	construct the Jacobian using previous array ("F" if filling by columns, "C", if filling by rows, $ROOTSYS/math/matrix/src/TMatrixT.cxx)
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG) << "Jacobian array converted to Jacobian matrix" << std::endl;

//	cluster covariance matrix by rows
	double totPFO_cov_matrix_by_rows[rows*rows] =
			{
				SigmaPx2		,	SigmaPxPy		,	SigmaPxPz		,	SigmaPxE	,
				SigmaPxPy		,	SigmaPy2		,	SigmaPyPz		,	SigmaPyE	,
				SigmaPxPz		,	SigmaPyPz		,	SigmaPz2		,	SigmaPzE	,
				SigmaPxE		,	SigmaPyE		,	SigmaPzE		,	SigmaE2
			};
	streamlog_out(DEBUG) << "cluster covariance matrix array formed by rows" << std::endl;

	TMatrixD covMatrix(rows,rows, totPFO_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG) << "cluster covariance matrix array converted to cluster covariance matrix" << std::endl;

	covMatrixMomenta.Mult( TMatrixD( jacobian ,
					TMatrixD::kTransposeMult ,
					covMatrix) ,
					jacobian
					);
	streamlog_out(DEBUG) << "cluster covariance matrix array in cartesian coordinate system converted to cluster covariance matrix array in spherical (polar) coordinate system" << std::endl;

	covP.push_back( covMatrixMomenta(0,0) ); // Theta-Theta
	covP.push_back( covMatrixMomenta(1,0) ); // Theta-Phi
	covP.push_back( covMatrixMomenta(1,1) ); // Phi-Phi
	covP.push_back( covMatrixMomenta(2,0) ); // Theta-P
	covP.push_back( covMatrixMomenta(2,1) ); // Phi-P
	covP.push_back( covMatrixMomenta(2,2) ); // P-P
	covP.push_back( covMatrixMomenta(3,0) ); // Theta-E
	covP.push_back( covMatrixMomenta(3,1) ); // Phi-E
	covP.push_back( covMatrixMomenta(3,2) ); // P-E
	covP.push_back( covMatrixMomenta(3,3) ); // E-E
	streamlog_out(DEBUG) << "FourMomentumCovarianceMatrix Filled succesfully" << std::endl;

	return covP;

}

TLorentzVector AddFourMomentumCovMatAllPFOs::getLinkedMCP( EVENT::LCEvent *pLCEvent, EVENT::ReconstructedParticle* inputPFO , int nTrackspfo , int nClusterspfo )
{
	LCRelationNavigator navClusterMCTruth(pLCEvent->getCollection(m_ClusterMCTruthLinkCollection));
	LCRelationNavigator navMCTruthCluster(pLCEvent->getCollection(m_MCTruthClusterLinkCollection));
	LCRelationNavigator navTrackMCTruth(pLCEvent->getCollection(m_TrackMCTruthLinkCollection));
	LCRelationNavigator navMCTruthTrack(pLCEvent->getCollection(m_MCTruthTrackLinkCollection));
	streamlog_out(DEBUG) << "PFO has " << nTrackspfo << " tracks and " << nClusterspfo << " clusters" << std::endl;

	int NeutralsPDGCode[14]{11,13,22,130,211,310,321,2112,2212,3112,3122,3222,3312,3322};
//	int ChargedPDGCode[14]{11,13,22,211,321,2212,3112,3122,3222,3312,3322};
	bool PFOlinkedtoMCP = false;
	TLorentzVector mcpFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	const EVENT::LCObjectVec& mcpvec = navClusterMCTruth.getRelatedToObjects(inputPFO->getClusters()[0]);
	const EVENT::FloatVec&  mcpweightvec = navClusterMCTruth.getRelatedToWeights(inputPFO->getClusters()[0]);
	MCParticle *linkedMCP;
	double maxweightPFOtoMCP = 0.0;
	int iPFOtoMCPmax = -1;
	int iMCPtoPFOmax = -1;
	streamlog_out(DEBUG) << "PFO is neutral (without track), pfoType: " << inputPFO->getType() << " , looking for linked " << mcpvec.size() << " MCPs" << std::endl;
	for ( unsigned int i_mcp = 0; i_mcp < mcpvec.size(); i_mcp++ )
	{
		double mcp_weight = mcpweightvec.at(i_mcp);
		MCParticle *testMCP = (MCParticle *) mcpvec.at(i_mcp);
		streamlog_out(DEBUG) << "checking linked MCP at " << i_mcp << " , MCP PDG = " << testMCP->getPDG() << " , link weight = " << mcp_weight << std::endl;
		if ( mcp_weight > maxweightPFOtoMCP && mcp_weight >= 0.9 )
		{
			maxweightPFOtoMCP = mcp_weight;
			iPFOtoMCPmax = i_mcp;
			streamlog_out(DEBUG) << "linkedMCP: " << i_mcp << " has PDG: " << testMCP->getPDG() << " and PFO to MCP Link has weight = " << mcp_weight << std::endl;
		}
	}
	if ( iPFOtoMCPmax != -1 )
	{
		h_NeutPFO_Weight->Fill( maxweightPFOtoMCP );
		linkedMCP = (MCParticle *) mcpvec.at(iPFOtoMCPmax);
		streamlog_out(DEBUG) << "Found linked MCP, MCP PDG: " << linkedMCP->getPDG() << " , link weight = " << maxweightPFOtoMCP << std::endl;
		Cluster *testCluster;
		const EVENT::LCObjectVec& clustervec = navMCTruthCluster.getRelatedToObjects(linkedMCP);
		const EVENT::FloatVec&  clusterweightvec = navMCTruthCluster.getRelatedToWeights(linkedMCP);
		double maxweightMCPtoPFO = 0.;
		for ( unsigned int i_cluster = 0; i_cluster < clustervec.size(); i_cluster++ )
		{
			double cluster_weight = clusterweightvec.at(i_cluster);
			testCluster = (Cluster *) clustervec.at(i_cluster);
			if ( cluster_weight > maxweightMCPtoPFO && cluster_weight >= 0.9 )
			{
				maxweightMCPtoPFO = cluster_weight;
				iMCPtoPFOmax = i_cluster;
			}
		}
		if ( iMCPtoPFOmax != -1 && testCluster == inputPFO->getClusters()[0] )
		{
			PFOlinkedtoMCP = true;
			h_NeutPFO_PDG->Fill( linkedMCP->getPDG() );
			if ( inputPFO->getType() != 22 ) streamlog_out(DEBUG) << "Initial PFO type: " << inputPFO->getType() << "	, linked MCP PDG(weight): " << linkedMCP->getPDG() << " (" << maxweightPFOtoMCP << ")	, linked-back PFO type(weight): " << inputPFO->getType() << " (" << maxweightMCPtoPFO << ")" << std::endl;
			bool KnownPFO = false;
			for ( int l = 0 ; l < 14 ; ++l)
			{
				if ( abs( linkedMCP->getPDG() ) == NeutralsPDGCode[ l ] )
				{
					h_NeutPFO_TYPE->Fill( l );
					KnownPFO = true;
					if ( inputPFO->getType() == 22 )
					{
						h_NeutPFO_IDasPhoton->Fill( l );
					}
					else
					{
						h_NeutPFO_IDasOther->Fill( l );
					}
				}
			}
			if ( !KnownPFO )
			{
				h_NeutPFO_TYPE->Fill( 14 );
				if ( inputPFO->getType() == 22 )
				{
					h_NeutPFO_IDasPhoton->Fill( 14 );
				}
				else
				{
					h_NeutPFO_IDasOther->Fill( 14 );
				}
			}
			mcpFourMomentum = TLorentzVector( linkedMCP->getMomentum()[0] , linkedMCP->getMomentum()[1] , linkedMCP->getMomentum()[2] , linkedMCP->getEnergy() );
			h_NHEnergy->Fill( mcpFourMomentum.E() , inputPFO->getEnergy() );
			h_NH_EclusterPlusMass_Emcp->Fill( mcpFourMomentum.E() , ( inputPFO->getClusters()[0] )->getEnergy() + inputPFO->getMass() );
		}
	}
	if ( PFOlinkedtoMCP )
	{
		return mcpFourMomentum;
	}
	else
	{
		return TLorentzVector( 0.0 , 0.0 , 0.0 , 0.0 );
	}

}

void AddFourMomentumCovMatAllPFOs::check(EVENT::LCEvent *pLCEvent)
{

	LCCollection *inputPfoCollection{};
	LCCollection *outputPfoCollection{};
	try
	{
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		outputPfoCollection = pLCEvent->getCollection(m_outputPfoCollection);
		int n_inputPFOs = inputPfoCollection->getNumberOfElements();
		int n_outputPFOs = outputPfoCollection->getNumberOfElements();
		streamlog_out(DEBUG) << " CHECK : processed events: " << m_nEvtSum << " (Number of inputPFOS: " << n_inputPFOs << " , Number of outputPFOs: " << n_outputPFOs <<")" << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input/Output collection not found in event " << m_nEvt << std::endl;
        }

}

void AddFourMomentumCovMatAllPFOs::end()
{

	m_pTFile->cd();
	m_pTTree->Write();
	m_Histograms->cd();
	h_clusterE_pfoE->Write();
	h_NeutPFO_PDG->Write();
	h_NeutPFO_TYPE->Write();
	h_NeutPFO_IDasPhoton->Write();
	h_NeutPFO_IDasOther->Write();
	h_NeutPFO_Weight->Write();
	h_NH_EclusterPlusMass_Emcp->Write();
	h_NHEnergy->Write();
	m_CovMatElements->cd();
	h_SigmaPx2->Write();
	h_SigmaPxPy->Write();
	h_SigmaPy2->Write();
	h_SigmaPxPz->Write();
	h_SigmaPyPz->Write();
	h_SigmaPz2->Write();
	h_SigmaPxE->Write();
	h_SigmaPyE->Write();
	h_SigmaPzE->Write();
	h_SigmaE2->Write();
	m_Photon->cd();
	h_ResidualEnergy_ph->Write();
	h_ResidualTheta_ph->Write();
	h_ResidualPhi_ph->Write();
	h_ErrorEnergy_ph->Write();
	h_ErrorTheta_ph->Write();
	h_ErrorPhi_ph->Write();
	h_NormalizedResidualEnergy_ph->Write();
	h_NormalizedResidualTheta_ph->Write();
	h_NormalizedResidualPhi_ph->Write();
	m_ErrorParameterization->cd();
	m_NeutralPFO->cd();
	h_ResidualEnergy_NH->Write();
	h_ResidualTheta_NH->Write();
	h_ResidualPhi_NH->Write();
	h_ErrorEnergy_NH->Write();
	h_ErrorTheta_NH->Write();
	h_ErrorPhi_NH->Write();
	h_NormalizedResidualEnergy_NH->Write();
	h_NormalizedResidualTheta_NH->Write();
	h_NormalizedResidualPhi_NH->Write();
	m_ErrorParameterization->cd();
	m_pTFile->Close();
	delete m_pTFile;

//	std::cout << " END : processed events: " << m_nEvtSum << std::endl;

}