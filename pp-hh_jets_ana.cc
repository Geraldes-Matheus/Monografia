#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"  
//#include "SlowJet.h"
#include <vector>
#include <cmath>

using namespace Pythia8;

int main(int argc, char* argv[]) {
    // Verifica se o arquivo LHE foi fornecido
    if (argc < 2) {
        std::cerr << "Uso: " << argv[0] << " arquivo.lhe" << std::endl;
        return 1;
    }

    // Inicializa o gerador Pythia
    Pythia pythia;

    // Configuração para ler o arquivo LHE
    pythia.readString("Beams:frameType = 4");  // Frame de LHE sem feixes
    pythia.readString("Beams:LHEF = " + std::string(argv[1]));  // Lê o arquivo LHE fornecido

    // Liga a geração de radiação inicial e final (ISR e FSR)
    pythia.readString("PartonLevel:ISR = on");
    pythia.readString("PartonLevel:FSR = on");
    pythia.readString("PartonLevel:MPI = on");   // Ativa as Interações Múltiplas de Partons (MPI)

    // Jatos
    pythia.readString("JetMatching:slowJetPower = -1"); // Jatos (anti-Kt)

    // Inicializa o Pythia
    if (!pythia.init()) {
        std::cerr << "Erro ao inicializar o Pythia!" << std::endl;
        return 1;
    }

    // Parâmetros para o detector de jatos SlowJet
    double radius = 0.4;
    double pTjetMin = 10.0;
    double etaMax = 5.0;
    int nSel = 2;

    // Configuração do detector de jatos SlowJet
    SlowJet slowJet(-1, radius, pTjetMin, etaMax, nSel, 1);

    // Criando histogramas com ROOT (TH1F)
    TH1F nJets("nJets", "Número de Jatos", 50, -0.5, 49.5);
    TH1F pTjets("pTjets", "pT dos Jatos", 100, 0., 500.);
    TH1F yJets("yJets", "Rapidez dos Jatos", 100, -5., 5.);
    TH1F phiJets("phiJets", "Ângulo Phi dos Jatos", 100, -M_PI, M_PI);
    TH1F distJets("distJets", "Distância entre Jatos", 100, 0., 10.);
    TH1F pTdiff("pTdiff", "Diferença de pT entre Jatos", 100, -100., 400.);
    TH1F invMassJets("invMassJets", "Massa Invariante entre Jatos", 100, 0., 500.);  // Novo histograma para a massa invariante
    TH1F numParticlesInJets("numParticlesInJets", "Número de Partículas por Jato", 50, -0.5, 49.5);  // Novo histograma para número de partículas

    // Configura o arquivo ROOT e a árvore
    TFile *file = TFile::Open("events4.root", "RECREATE");
    if (!file) {
        std::cerr << "Erro ao abrir o arquivo ROOT!" << std::endl;
        return 1;
    }

    TTree *tree = new TTree("Events", "Eventos gerados pelo Pythia");

    // Variáveis para armazenar características das partículas
    std::vector<float> pt, eta, phi, mass, charge, energy;
    std::vector<float> px, py, pz;  // Componentes do momento
    std::vector<int> pdgId, status;

    // Cria as branches para as variáveis
    tree->Branch("pt", &pt);
    tree->Branch("eta", &eta);
    tree->Branch("phi", &phi);
    tree->Branch("mass", &mass);
    tree->Branch("charge", &charge);
    tree->Branch("energy", &energy);  // Energia
    tree->Branch("px", &px);            // Momento x
    tree->Branch("py", &py);            // Momento y
    tree->Branch("pz", &pz);            // Momento z
    tree->Branch("pdgId", &pdgId);
    tree->Branch("status", &status);

    // Loop para eventos
    for (int iEvent = 0; iEvent < 10000; ++iEvent) {
        // Processa o próximo evento
        if (!pythia.next()) continue;

        // Limpa as variáveis para o novo evento
        pt.clear();
        eta.clear();
        phi.clear();
        mass.clear();
        charge.clear();
        energy.clear();
        px.clear();
        py.clear();
        pz.clear();
        pdgId.clear();
        status.clear();

        // Analisa as propriedades dos jatos com SlowJet
        slowJet.analyze(pythia.event);
        if (iEvent < 10) slowJet.list();  // Mostra a lista dos primeiros jatos

        // Preenche os histogramas para as distribuições dos jatos
        nJets.Fill(slowJet.sizeJet());
        for (int i = 0; i < slowJet.sizeJet(); ++i) {
            pTjets.Fill(slowJet.pT(i));
            yJets.Fill(slowJet.y(i));
            phiJets.Fill(slowJet.phi(i));

            // Conta o número de partículas no jato
            int numParticles = 0;
            for (int j = 0; j < pythia.event.size(); ++j) {
                Particle& particle = pythia.event[j];
                double dEta = std::abs(slowJet.y(i) - particle.eta());
                double dPhi = std::abs(slowJet.phi(i) - particle.phi());
                if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
                double dR = std::sqrt(dEta * dEta + dPhi * dPhi);
                if (dR < radius && particle.pT() > pTjetMin) {  // Verifica se a partícula está dentro do cone do jato
                    numParticles++;
                }
            }
            numParticlesInJets.Fill(numParticles);  // Preenche o histograma com o número de partículas no jato
        }

        // Preenche o histograma de distância entre jatos
        for (int i = 0; i < slowJet.sizeJet() - 1; ++i) {
            for (int j = i + 1; j < slowJet.sizeJet(); ++j) {
                double dY = slowJet.y(i) - slowJet.y(j);
                double dPhi = std::abs(slowJet.phi(i) - slowJet.phi(j));
                if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
                double dR = std::sqrt(dY * dY + dPhi * dPhi);
                distJets.Fill(dR);

                // Calculando a massa invariante entre o par de jatos (i, j) usando a fórmula corrigida
                double pT1 = slowJet.pT(i);
                double pT2 = slowJet.pT(j);
                double eta1 = slowJet.y(i);  // Rapidez do primeiro jato
                double eta2 = slowJet.y(j);  // Rapidez do segundo jato
                double phi1 = slowJet.phi(i);  // Ângulo phi do primeiro jato
                double phi2 = slowJet.phi(j);  // Ângulo phi do segundo jato

                // Usando a fórmula corrigida para a massa invariante
                double invMass2 = 2 * pT1 * pT2 * (std::cosh(eta1 - eta2) - std::cos(phi1 - phi2));

                if (invMass2 > 0) {
                    invMassJets.Fill(std::sqrt(invMass2));  // Massa invariante real (se for positiva)
                }
            }
        }

        // Preenche o histograma da diferença de pT entre jatos
        for (int i = 1; i < slowJet.sizeJet(); ++i) {
            pTdiff.Fill(slowJet.pT(i-1) - slowJet.pT(i));
        }

        // Preenche a árvore com as partículas do evento
        for (int i = 0; i < pythia.event.size(); ++i) {
            Particle& particle = pythia.event[i];
            pt.push_back(particle.pT());
            eta.push_back(particle.eta());
            phi.push_back(particle.phi());
            mass.push_back(particle.m());
            charge.push_back(particle.charge());
            energy.push_back(particle.e());
            px.push_back(particle.px());
            py.push_back(particle.py());
            pz.push_back(particle.pz());
            pdgId.push_back(particle.id());
            status.push_back(particle.status());
        }

        // Preenche a árvore com o evento atual
        tree->Fill();
    }

    // Estatísticas de fim da execução
    pythia.stat();

    // Salva os histogramas no arquivo ROOT
    file->cd();  // Certifique-se de que estamos no diretório correto dentro do arquivo ROOT
    nJets.Write();  // Salva os histogramas no arquivo
    pTjets.Write();
    yJets.Write();
    phiJets.Write();
    distJets.Write();
    pTdiff.Write();
    invMassJets.Write();  // Salva o histograma de massa invariante
    numParticlesInJets.Write();  // Salva o histograma de número de partículas nos jatos

    // Salva a árvore de eventos
    tree->Write();
    delete file;  // Fecha o arquivo ROOT

    // Done.
    return 0;
}
