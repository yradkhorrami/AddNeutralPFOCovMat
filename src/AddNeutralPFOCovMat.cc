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

AddNeutralPFOCovMat aAddNeutralPFOCovMat;

AddNeutralPFOCovMat::AddNeutralPFOCovMat() :

Processor("AddNeutralPFOCovMat"),
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

	registerProcessorParameter(	"useTrueJacobian",
					"true: Use (mathematically) true Jacobian for the option E_cluster = |p|, false: for the option E_cluster = |p|, Use the same jacobian as the option E_cluster = E_kinetic",
					m_useTrueJacobian,
					bool(false)
				);

	registerProcessorParameter(	"MinWeightClusterMCTruthLink" ,
					"Minimum acceptable weight for Cluster -> MCParticle Link"  ,
					m_MinWeightClusterMCTruthLink ,
					float(0.9f)
				);

	registerProcessorParameter(	"MinWeightMCTruthClusterLink" ,
					"Minimum acceptable weight for MCParticle -> Cluster Link"  ,
					m_MinWeightMCTruthClusterLink ,
					float(0.9f)
				);

	registerProcessorParameter(	"fillRootTree",
					"whether store output comparison in RootTree or not",
					m_fillRootTree,
					bool(false)
				);

	registerProcessorParameter(	"RootFile",
	                                "Name of the output root file",
					m_rootFile,
					std::string("FourMomentumCovMatAllPFOs.root")
				);

}

void AddNeutralPFOCovMat::init()
{

	streamlog_out(MESSAGE) << "   init called  " << std::endl;
	printParameters();
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

void AddNeutralPFOCovMat::Clear()
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

void AddNeutralPFOCovMat::processRunHeader()
{

	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;

}

void AddNeutralPFOCovMat::processEvent( EVENT::LCEvent *pLCEvent )
{

	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	++m_nEvtSum;

	LCCollection *inputPfoCollection{};
	LCCollectionVec *outputPfoCollection{};
	int n_PFO = -1;
	this->Clear();
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;

	try
	{
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		outputPfoCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
		n_PFO = inputPfoCollection->getNumberOfElements();
		if ( n_PFO == -1 ) streamlog_out(DEBUG7) << "	Input PFO collection (" << m_inputPfoCollection << ") has no element (PFO) " << std::endl;
		streamlog_out(DEBUG7) << "	Total Number of PFOs: " << n_PFO << std::endl;
		for (int i_pfo = 0; i_pfo < n_PFO ; ++i_pfo)
		{
			streamlog_out(DEBUG6) << "	-------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG6) << "	Processing PFO at index " << i_pfo << std::endl;
			streamlog_out(DEBUG6) << "" << std::endl;
			ReconstructedParticle* inputPFO = dynamic_cast<ReconstructedParticle*>( inputPfoCollection->getElementAt( i_pfo ) );
			ReconstructedParticleImpl* outputPFO = new ReconstructedParticleImpl;
			int linkedMCP_PDGCode = 0;
			float weightClusterMCP = 0.0;
			float weightMCPCluster = 0.0;
			bool m_updatePFO = true;
			float pfoMass = inputPFO->getMass();
			TVector3 clusterPosition( 0.0 , 0.0 , 0.0 );
			TLorentzVector mcpFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			TLorentzVector pfoFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
			double outputPFOMomentum[3]{0., 0., 0.};
			std::vector<float> outputCovMatrix( 10 , 0.0 );
			std::vector<float> PFOResidual( 3 , 0.0 );
			std::vector<float> PFOCovMatPolar( 10 , 0.0 );
			std::vector<float> PFOCoordinateError( 6 , 0.0 );
			if ( ( inputPFO->getTracks() ).size() !=0 )
			{
				streamlog_out(DEBUG) << "	PFO has one (or more) track(s), Track parameters are used for PFO CovMat. nothing to do/update!" << std::endl;
				m_updatePFO = false;
				outputPFO->setType(inputPFO->getType());
				outputPFO->setMomentum( inputPFO->getMomentum() );
				outputPFO->setEnergy( inputPFO->getEnergy() );
				outputPFO->setMass( inputPFO->getMass() );
				outputPFO->setCovMatrix(inputPFO->getCovMatrix());
				outputPFO->setCharge(inputPFO->getCharge());
				outputPFO->setReferencePoint(inputPFO->getReferencePoint());
				for (unsigned int j=0; j<inputPFO->getParticleIDs().size(); ++j)
				{
					ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>(inputPFO->getParticleIDs()[j]);
				        ParticleIDImpl* outPID = new ParticleIDImpl;
				        outPID->setType(inPID->getType());
				        outPID->setPDG(inPID->getPDG());
				        outPID->setLikelihood(inPID->getLikelihood());
				        outPID->setAlgorithmType(inPID->getAlgorithmType()) ;
				        for (unsigned int k=0; k<inPID->getParameters().size()  ; ++k) outPID->addParameter(inPID->getParameters()[k]) ;
				        outputPFO->addParticleID(outPID);
				}
				outputPFO->setParticleIDUsed(inputPFO->getParticleIDUsed());
				outputPFO->setGoodnessOfPID(inputPFO->getGoodnessOfPID());
				for (unsigned int j=0; j<inputPFO->getParticles().size(); ++j)
				{
					outputPFO->addParticle(inputPFO->getParticles()[j]);
				}
				for (unsigned int j=0; j<inputPFO->getClusters().size(); ++j)
				{
					outputPFO->addCluster(inputPFO->getClusters()[j]);
				}
				for (unsigned int j=0; j<inputPFO->getTracks().size(); ++j)
				{
					outputPFO->addTrack(inputPFO->getTracks()[j]);
				}
				outputPFO->setStartVertex(inputPFO->getStartVertex());
			}
			else
			{
				mcpFourMomentum = this->getLinkedMCP( pLCEvent , inputPFO, linkedMCP_PDGCode , weightClusterMCP , weightMCPCluster );
				streamlog_out(DEBUG2) << "	PDG code of linked MCParticle is: " << linkedMCP_PDGCode << "( " << weightClusterMCP << " , " << weightMCPCluster << " )" << std::endl;
				if ( !m_AssumeNeutralPFOMassive ) pfoMass = 0.0;
				float clusterEnergy	= ( inputPFO->getClusters()[0] )->getEnergy();
				float clusterX		= ( inputPFO->getClusters()[0] )->getPosition()[0];
				float clusterY		= ( inputPFO->getClusters()[0] )->getPosition()[1];
				float clusterZ		= ( inputPFO->getClusters()[0] )->getPosition()[2];
				clusterPosition	= TVector3( clusterX , clusterY , clusterZ );
				float clusterDistance	= sqrt( pow( clusterX , 2 ) + pow( clusterY , 2 ) + pow( clusterZ , 2 ) );
				float pfoMomentumMag	= 0;
				float pfoEnergy	= inputPFO->getEnergy();
				float pfoE;
				pfoMomentumMag = ( m_isClusterEnergyKinEnergy ? sqrt( pow( pfoEnergy , 2 ) + 2 * pfoMass * pfoEnergy ) : pfoEnergy );
				float pfoPx;
				float pfoPy;
				float pfoPz;
				if ( m_updatePFO4Momentum )
				{
					pfoPx	= pfoMomentumMag * clusterX / clusterDistance;
					pfoPy	= pfoMomentumMag * clusterY / clusterDistance;
					pfoPz	= pfoMomentumMag * clusterZ / clusterDistance;
					pfoE = ( m_isClusterEnergyKinEnergy ? pfoEnergy + pfoMass : sqrt( pow( pfoMomentumMag , 2 ) + pow( pfoMass , 2 ) ) );
				}
				else
				{
					pfoPx	= inputPFO->getMomentum()[ 0 ];
					pfoPy	= inputPFO->getMomentum()[ 1 ];
					pfoPz	= inputPFO->getMomentum()[ 2 ];
					pfoE	= inputPFO->getEnergy();
				}
				std::vector<float> clusterPositionError = ( inputPFO->getClusters()[0] )->getPositionError();
				float clusterEnergyError = ( inputPFO->getClusters()[0] )->getEnergyError();
				streamlog_out(DEBUG2) << "	Cluster Energy / PFO Energy = " << clusterEnergy << " / " << pfoE << std::endl;
				h_clusterE_pfoE->Fill( clusterEnergy , pfoE );
				TVector3 pfoMomentum( pfoPx , pfoPy , pfoPz );
				pfoFourMomentum	= TLorentzVector( pfoMomentum , pfoE );
				m_RecoEnergy.push_back( pfoFourMomentum.E() );
				m_RecoTheta.push_back( pfoFourMomentum.Theta() );
				m_RecoPhi.push_back( pfoFourMomentum.Phi() );
				outputCovMatrix	= this->UpdateNeutralPFOCovMat( clusterPosition , pfoEnergy , pfoMass , clusterPositionError , clusterEnergyError );
				outputPFOMomentum[ 0 ] = pfoFourMomentum.Px();
				outputPFOMomentum[ 1 ] = pfoFourMomentum.Py();
				outputPFOMomentum[ 2 ] = pfoFourMomentum.Pz();
				pfoE = pfoFourMomentum.E();

				outputPFO->setType(inputPFO->getType());
				outputPFO->setMomentum( outputPFOMomentum );
				outputPFO->setEnergy( pfoE );
				outputPFO->setMass( inputPFO->getMass() );
				outputPFO->setCovMatrix( outputCovMatrix );
				outputPFO->setCharge(inputPFO->getCharge());
				outputPFO->setReferencePoint(inputPFO->getReferencePoint());
				for (unsigned int j=0; j<inputPFO->getParticleIDs().size(); ++j)
				{
					ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>(inputPFO->getParticleIDs()[j]);
				        ParticleIDImpl* outPID = new ParticleIDImpl;
				        outPID->setType(inPID->getType());
				        outPID->setPDG(inPID->getPDG());
				        outPID->setLikelihood(inPID->getLikelihood());
				        outPID->setAlgorithmType(inPID->getAlgorithmType()) ;
				        for (unsigned int k=0; k<inPID->getParameters().size()  ; ++k) outPID->addParameter(inPID->getParameters()[k]) ;
				        outputPFO->addParticleID(outPID);
				}
				outputPFO->setParticleIDUsed(inputPFO->getParticleIDUsed());
				outputPFO->setGoodnessOfPID(inputPFO->getGoodnessOfPID());
				for (unsigned int j=0; j<inputPFO->getParticles().size(); ++j)
				{
					outputPFO->addParticle(inputPFO->getParticles()[j]);
				}
				for (unsigned int j=0; j<inputPFO->getClusters().size(); ++j)
				{
					outputPFO->addCluster(inputPFO->getClusters()[j]);
				}
				for (unsigned int j=0; j<inputPFO->getTracks().size(); ++j)
				{
					outputPFO->addTrack(inputPFO->getTracks()[j]);
				}
				outputPFO->setStartVertex(inputPFO->getStartVertex());
			}
			outputPfoCollection->addElement( outputPFO );
		}
		pLCEvent->addCollection( outputPfoCollection , m_outputPfoCollection );
		m_pTTree->Fill();
	}
        catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input collection not found in event " << m_nEvt << std::endl;
        }

}

std::vector<float> AddNeutralPFOCovMat::UpdateNeutralPFOCovMat( TVector3 clusterPosition , float pfoEc , float pfoMass , std::vector<float> clusterPositionError , float clusterEnergyError )
{

//	Obtain covariance matrix on (px,py,pz,E) from the
//	covariance matrix on cluster parameters (px,py,pz,|p|=Ec).
//	=> E^2 = Ec^2 + m^2	;	|p| = Ec
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
//			Dpx/DEc			Dpy/DEc			Dpz/DEc			DE/DEc
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

	float pfoP = ( m_isClusterEnergyKinEnergy ? sqrt( pow( pfoEc , 2 ) + 2 * pfoMass * pfoEc ) : pfoEc );
	float pfoE = ( m_isClusterEnergyKinEnergy ? sqrt( pow( pfoP , 2 ) + pow( pfoMass , 2 ) ) : sqrt( pow( pfoP , 2 ) + pow( pfoMass , 2 ) ) );
	float derivative_coeff	= ( ( m_useTrueJacobian && !m_isClusterEnergyKinEnergy ) ? pfoP / pfoE : 1.0 );

	streamlog_out(DEBUG0) << "	Cluster information obtained" << std::endl;

//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
	{
		pfoP * ( pfoR2 - pfoX2 ) / pfoR3			,	-pfoP * pfoX * pfoY / pfoR3				,	-pfoP * pfoX * pfoZ / pfoR3				,	0			,
		-pfoP * pfoY * pfoX / pfoR3				,	pfoP * ( pfoR2 - pfoY2 ) / pfoR3			,	-pfoP * pfoY * pfoZ / pfoR3				,	0			,
		-pfoP * pfoZ * pfoX / pfoR3				,	-pfoP * pfoZ * pfoY / pfoR3				,	pfoP * ( pfoR2 - pfoZ2 ) / pfoR3			,	0			,
		derivative_coeff * pfoE * pfoX / ( pfoP * pfoR )	,	derivative_coeff * pfoE * pfoY / ( pfoP * pfoR )	,	derivative_coeff * pfoE * pfoZ / ( pfoP * pfoR )	,	derivative_coeff
	};

//	construct the Jacobian using previous array ("F" if filling by columns, "C" if filling by rows, $ROOTSYS/math/matrix/src/TMatrixT.cxx)
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG0) << "	Jacobian array converted to Jacobian matrix" << std::endl;

//	cluster covariance matrix by rows
	double cluster_cov_matrix_by_rows[rows*rows] =
			{
				SigmaX2		,	SigmaXY		,	SigmaXZ		,	0	,
				SigmaXY		,	SigmaY2		,	SigmaYZ		,	0	,
				SigmaXZ		,	SigmaYZ		,	SigmaZ2		,	0	,
				0		,	0		,	0		,	SigmaE2
			};

	TMatrixD covMatrix_cluster(rows,rows, cluster_cov_matrix_by_rows, "C");
	streamlog_out(DEBUG0) << "	Cluster covariance matrix array converted to cluster covariance matrix" << std::endl;

	covMatrixMomenta.Mult( TMatrixD( jacobian ,
					TMatrixD::kTransposeMult ,
					covMatrix_cluster) ,
					jacobian
					);

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
	streamlog_out(DEBUG0) << "	FourMomentumCovarianceMatrix Filled succesfully" << std::endl;

	return covP;

}

std::vector<float> AddNeutralPFOCovMat::getClusterDirectionError( TVector3 clusterPosition , std::vector<float> clusterPositionError )
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

std::vector<float> AddNeutralPFOCovMat::getPFOResidual( TLorentzVector pfoFourMomentum , TLorentzVector mcpFourMomentum )
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

std::vector<float> AddNeutralPFOCovMat::getPFOCovMatPolarCoordinate( TLorentzVector pfoFourMomentum , std::vector<float> pfoCovMat )
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

TLorentzVector AddNeutralPFOCovMat::getLinkedMCP( EVENT::LCEvent *pLCEvent, EVENT::ReconstructedParticle* inputPFO , int &linkedMCP_PDGCode , float &weightClusterMCP , float &weightMCPCluster )
{
	LCRelationNavigator navClusterMCTruth(pLCEvent->getCollection(m_ClusterMCTruthLinkCollection));
	LCRelationNavigator navMCTruthCluster(pLCEvent->getCollection(m_MCTruthClusterLinkCollection));

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
	streamlog_out(DEBUG0) << "	PFO is neutral (without track), pfoType: " << inputPFO->getType() << " , looking for linked " << mcpvec.size() << " MCPs" << std::endl;
	for ( unsigned int i_mcp = 0; i_mcp < mcpvec.size(); i_mcp++ )
	{
		double mcp_weight = mcpweightvec.at(i_mcp);
		MCParticle *testMCP = (MCParticle *) mcpvec.at(i_mcp);
		streamlog_out(DEBUG0) << "	Checking MCP[ " << i_mcp << " ] , MCP PDG = " << testMCP->getPDG() << " , link weight = " << mcp_weight << std::endl;
		if ( mcp_weight > maxweightPFOtoMCP )// && mcp_weight >= 0.9 )
		{
			maxweightPFOtoMCP = mcp_weight;
			iPFOtoMCPmax = i_mcp;
		}
	}
	weightClusterMCP = maxweightPFOtoMCP;
	if ( iPFOtoMCPmax != -1 )
	{
		h_NeutPFO_Weight->Fill( maxweightPFOtoMCP );
		linkedMCP = (MCParticle *) mcpvec.at( iPFOtoMCPmax );
		Cluster *linkedCluster;
		const EVENT::LCObjectVec& clustervec = navMCTruthCluster.getRelatedToObjects(linkedMCP);
		const EVENT::FloatVec&  clusterweightvec = navMCTruthCluster.getRelatedToWeights(linkedMCP);
		streamlog_out(DEBUG0) << "	Found linked MCP, MCP PDG: " << linkedMCP->getPDG() << " , link weight = " << maxweightPFOtoMCP << " , looking for " << clustervec.size() << " Cluster(s) linked back to this MCParticle" << std::endl;
		double maxweightMCPtoPFO = 0.;
		for ( unsigned int i_cluster = 0; i_cluster < clustervec.size(); i_cluster++ )
		{
			double cluster_weight = clusterweightvec.at(i_cluster);
			streamlog_out(DEBUG0) << "	Checking cluster[ " << i_cluster << " ] , link weight = " << cluster_weight << std::endl;
			if ( cluster_weight > maxweightMCPtoPFO )// && cluster_weight >= 0.9 )
			{
				maxweightMCPtoPFO = cluster_weight;
				iMCPtoPFOmax = i_cluster;
			}
		}
		weightMCPCluster = maxweightMCPtoPFO;
		if ( iMCPtoPFOmax != -1 )
		{
			linkedCluster = (Cluster *) clustervec.at( iMCPtoPFOmax );
			if ( linkedCluster == inputPFO->getClusters()[0] )
			{
				streamlog_out(DEBUG1) << "	Found a MCParticle (PDGCode: " << linkedMCP->getPDG() << ") 	linked to cluster of PFO (TYPE: " << inputPFO->getType() << ")" << std::endl;
				streamlog_out(DEBUG1) << "	Cluster_MCParticle link weight = " << maxweightPFOtoMCP << " 	; 	MCParticle_Cluster link weight = " << maxweightMCPtoPFO << std::endl;
				PFOlinkedtoMCP = true;
				linkedMCP_PDGCode = linkedMCP->getPDG();
				h_NeutPFO_PDG->Fill( linkedMCP->getPDG() );
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

void AddNeutralPFOCovMat::check(EVENT::LCEvent *pLCEvent)
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

void AddNeutralPFOCovMat::end()
{

	if ( m_fillRootTree )
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
		m_pTFile->Close();
		delete m_pTFile;
	}

//	std::cout << " END : processed events: " << m_nEvtSum << std::endl;

}
