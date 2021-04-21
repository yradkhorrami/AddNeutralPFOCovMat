#ifndef AddNeutralPFOCovMat_h
#define AddNeutralPFOCovMat_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "UTIL/LCRelationNavigator.h"
#include <EVENT/MCParticle.h>
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <vector>
#include "TMatrixD.h"
#include <TFile.h>
#include <TTree.h>

class TFile;
class TDirectory;
class TH1F;
class TH1I;
class TH2I;
class TH2F;
class TTree;

using namespace lcio ;
using namespace marlin ;

class AddNeutralPFOCovMat : public Processor
{
	public:

		virtual Processor*  newProcessor()
		{
			return new AddNeutralPFOCovMat;
		}
		AddNeutralPFOCovMat();
		virtual ~AddNeutralPFOCovMat() = default;
		AddNeutralPFOCovMat(const AddNeutralPFOCovMat&) = delete;
		AddNeutralPFOCovMat& operator=(const AddNeutralPFOCovMat&) = delete;
		virtual void init();
		virtual void Clear();
		virtual void processRunHeader();
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		std::vector<float> UpdateNeutralPFOCovMat( TVector3 clusterPosition , float pfoEc , float pfoMass , std::vector<float> clusterPositionError , float clusterEnergyError );
		std::vector<float> getClusterDirectionError( TVector3 clusterPosition , std::vector<float> clusterPositionError );
		std::vector<float> getPFOResidual( TLorentzVector pfoFourMomentum , TLorentzVector mcpFourMomentum );
		std::vector<float> getPFOCovMatPolarCoordinate( TLorentzVector pfoFourMomentum , std::vector<float> pfoCovMat );
		TLorentzVector getLinkedMCP( EVENT::LCEvent *pLCEvent , EVENT::ReconstructedParticle* inputPFO , int &linkedMCP_PDGCode , float &weightClusterMCP , float &weightMCPCluster );
		virtual void check( EVENT::LCEvent *pLCEvent );
		virtual void end();

	private:
		std::string				m_inputPfoCollection{};
		std::string				m_ClusterMCTruthLinkCollection{};
		std::string				m_MCTruthClusterLinkCollection{};
		std::string				m_outputPfoCollection{};
		std::string				m_rootFile{};

		bool					m_AssumeNeutralPFOMassive = true;
		bool					m_isClusterEnergyKinEnergy = false;
		bool					m_updatePFO4Momentum = false;
		bool					m_useTrueJacobian = false;
		bool					m_fillRootTree = false;
		float					m_MinWeightClusterMCTruthLink;
		float					m_MinWeightMCTruthClusterLink;
		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;
		TFile					*m_pTFile;
	        TTree					*m_pTTree;
		std::vector<int>			m_foundLinkedMCP{};
		std::vector<float>			m_mcEnergy{};
		std::vector<float>			m_mcTheta{};
		std::vector<float>			m_mcPhi{};
		std::vector<float>			m_RecoEnergy{};
		std::vector<float>			m_RecoTheta{};
		std::vector<float>			m_RecoPhi{};
		std::vector<float>			m_ResidualEnergy{};
		std::vector<float>			m_ResidualTheta{};
		std::vector<float>			m_ResidualPhi{};
		std::vector<float>			m_ErrorEnergy{};
		std::vector<float>			m_ErrorTheta{};
		std::vector<float>			m_ErrorPhi{};
		std::vector<float>			m_NormalizedResidualEnergy{};
		std::vector<float>			m_NormalizedResidualTheta{};
		std::vector<float>			m_NormalizedResidualPhi{};
		std::vector<int>			m_foundLinkedMCP_Ph{};
		std::vector<float>			m_mcEnergy_Ph{};
		std::vector<float>			m_mcTheta_Ph{};
		std::vector<float>			m_mcPhi_Ph{};
		std::vector<float>			m_RecoEnergy_Ph{};
		std::vector<float>			m_RecoTheta_Ph{};
		std::vector<float>			m_RecoPhi_Ph{};
		std::vector<float>			m_ResidualEnergy_Ph{};
		std::vector<float>			m_ResidualTheta_Ph{};
		std::vector<float>			m_ResidualPhi_Ph{};
		std::vector<float>			m_ErrorEnergy_Ph{};
		std::vector<float>			m_ErrorTheta_Ph{};
		std::vector<float>			m_ErrorPhi_Ph{};
		std::vector<float>			m_NormalizedResidualEnergy_Ph{};
		std::vector<float>			m_NormalizedResidualTheta_Ph{};
		std::vector<float>			m_NormalizedResidualPhi_Ph{};
		std::vector<int>			m_foundLinkedMCP_NH{};
		std::vector<float>			m_mcEnergy_NH{};
		std::vector<float>			m_mcTheta_NH{};
		std::vector<float>			m_mcPhi_NH{};
		std::vector<float>			m_RecoEnergy_NH{};
		std::vector<float>			m_RecoTheta_NH{};
		std::vector<float>			m_RecoPhi_NH{};
		std::vector<float>			m_ResidualEnergy_NH{};
		std::vector<float>			m_ResidualTheta_NH{};
		std::vector<float>			m_ResidualPhi_NH{};
		std::vector<float>			m_ErrorEnergy_NH{};
		std::vector<float>			m_ErrorTheta_NH{};
		std::vector<float>			m_ErrorPhi_NH{};
		std::vector<float>			m_NormalizedResidualEnergy_NH{};
		std::vector<float>			m_NormalizedResidualTheta_NH{};
		std::vector<float>			m_NormalizedResidualPhi_NH{};
	        TDirectory				*m_Histograms;
	        TDirectory				*m_CovMatElements;
	        TDirectory				*m_NeutralPFOswithoutTrak;
	        TDirectory				*m_Photon;
	        TDirectory				*m_NeutralPFO;
		TH2F					*h_clusterE_pfoE{};
		TH2F					*h_SigmaPx2{};
		TH2F					*h_SigmaPxPy{};
		TH2F					*h_SigmaPy2{};
		TH2F					*h_SigmaPxPz{};
		TH2F					*h_SigmaPyPz{};
		TH2F					*h_SigmaPz2{};
		TH2F					*h_SigmaPxE{};
		TH2F					*h_SigmaPyE{};
		TH2F					*h_SigmaPzE{};
		TH2F					*h_SigmaE2{};
		TH1I					*h_NeutPFO_PDG{};
		TH1I					*h_NeutPFO_TYPE{};
		TH1I					*h_NeutPFO_IDasPhoton{};
		TH1I					*h_NeutPFO_IDasOther{};
		TH1F					*h_NeutPFO_Weight{};
		TH1F					*h_ResidualEnergy_ph{};
		TH1F					*h_ResidualTheta_ph{};
		TH1F					*h_ResidualPhi_ph{};
		TH1F					*h_ErrorEnergy_ph{};
		TH1F					*h_ErrorTheta_ph{};
		TH1F					*h_ErrorPhi_ph{};
		TH1F					*h_NormalizedResidualEnergy_ph{};
		TH1F					*h_NormalizedResidualTheta_ph{};
		TH1F					*h_NormalizedResidualPhi_ph{};
		TH1F					*h_ResidualEnergy_NH{};
		TH1F					*h_ResidualTheta_NH{};
		TH1F					*h_ResidualPhi_NH{};
		TH1F					*h_ErrorEnergy_NH{};
		TH1F					*h_ErrorTheta_NH{};
		TH1F					*h_ErrorPhi_NH{};
		TH1F					*h_NormalizedResidualEnergy_NH{};
		TH1F					*h_NormalizedResidualTheta_NH{};
		TH1F					*h_NormalizedResidualPhi_NH{};
		TH2F					*h_NH_EclusterPlusMass_Emcp{};
		TH2F					*h_NHEnergy{};

};
#endif
