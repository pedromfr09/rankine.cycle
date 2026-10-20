"""
Simulador de ciclo de Rankine com irreversibilidades (água como fluido de trabalho).

O usuário informa as pressões, a temperatura de entrada na turbina e as
eficiências isentrópicas da bomba e da turbina. O programa calcula os quatro
estados do ciclo com o CoolProp, mostra o desempenho e gera o diagrama T-s.

Hipóteses:
  - Estado 1: líquido saturado na pressão do condensador.
  - Bomba: líquido incompressível (v constante, avaliado no estado 1).
  - Caldeira e condensador: isobáricos, sem perda de carga.
  - Turbina e bomba: adiabáticas, com eficiência isentrópica.
  - Variações de energia cinética e potencial desprezadas.
"""

import matplotlib.pyplot as plt
import numpy as np
from CoolProp.CoolProp import PropsSI


def coletar_medicao(nome_parametro, unidade, limite_inferior=0.0, limite_superior=None):
    """Pede um número ao usuário até receber um valor válido.

    Aceita ponto ou vírgula como separador decimal. O valor deve ser maior
    que limite_inferior e, se informado, menor ou igual a limite_superior.
    """
    while True:
        dado_bruto = input(f"medição no equipamento - {nome_parametro} ({unidade}): ")

        try:
            valor = float(dado_bruto.replace(",", "."))
        except ValueError:
            print("ERRO: use apenas números (ex: 45.5 ou 45,5).")
            continue

        if not np.isfinite(valor):
            print("ERRO: informe um número finito.")
        elif valor <= limite_inferior:
            print(f"ERRO: o valor deve ser maior que {limite_inferior}.")
        elif limite_superior is not None and valor > limite_superior:
            print(f"ERRO: o valor deve ser menor ou igual a {limite_superior}.")
        else:
            return valor


def linha_isobarica(fluido, P, h_inicial, h_final, n=100):
    """Pontos (s em kJ/kg.K, T em °C) de uma isobárica entre duas entalpias.

    Usar a entalpia como variável evita problemas no patamar de mudança
    de fase, onde T é constante e só s e h variam.
    """
    h = np.linspace(h_inicial, h_final, n)
    s = [PropsSI('S', 'P', P, 'H', hi, fluido) / 1000 for hi in h]
    T = [PropsSI('T', 'P', P, 'H', hi, fluido) - 273.15 for hi in h]
    return s, T


def simular_rankine():
    print("\n" + "=" * 55)
    print("      simulador de ciclo de rankine")
    print("=" * 55)

    # coleta de dados
    P_cond_bar = coletar_medicao("pressão do Condensador", "bar")
    P_cald_bar = coletar_medicao("pressão da Caldeira", "bar")
    T_turb_C = coletar_medicao("temperatura de entrada na Turbina", "°C")
    eta_bomba = coletar_medicao("eficiência Isentrópica da Bomba", "entre 0 e 1, ex: 0.85",
                                limite_superior=1.0)
    eta_turb = coletar_medicao("eficiência Isentrópica da Turbina", "entre 0 e 1, ex: 0.90",
                               limite_superior=1.0)

    P_baixa = P_cond_bar * 1e5   # bar -> Pa
    P_alta = P_cald_bar * 1e5
    T_3 = T_turb_C + 273.15      # °C -> K
    fluido = 'Water'

    # validações físicas antes de calcular
    P_crit = PropsSI('Pcrit', fluido)

    if P_alta <= P_baixa:
        print("\nERRO: a pressão da caldeira deve ser maior que a do condensador.")
        return

    if P_baixa >= P_crit:
        print(f"\nERRO: a pressão do condensador deve ser menor que a crítica "
              f"({P_crit / 1e5:.1f} bar).")
        return

    if P_alta < P_crit:
        T_sat_alta = PropsSI('T', 'P', P_alta, 'Q', 1, fluido)
        if T_3 <= T_sat_alta:
            print(f"\nERRO: {T_turb_C:.1f} °C não é vapor superaquecido a {P_cald_bar:.1f} bar "
                  f"(T_sat = {T_sat_alta - 273.15:.1f} °C). "
                  f"Informe uma temperatura maior que a de saturação.")
            return

    print("\ncalculando propriedades termofísicas via CoolProp")

    try:
        # estado 1: entrada da bomba (líquido saturado no condensador)
        h1 = PropsSI('H', 'P', P_baixa, 'Q', 0, fluido)
        s1 = PropsSI('S', 'P', P_baixa, 'Q', 0, fluido)
        T1 = PropsSI('T', 'P', P_baixa, 'Q', 0, fluido)
        v1 = 1 / PropsSI('D', 'P', P_baixa, 'Q', 0, fluido)  # D é densidade (kg/m³)

        # estado 2: saída da bomba (líquido sub-resfriado)
        # líquido incompressível: w = v * ΔP
        w_bomba_ideal = v1 * (P_alta - P_baixa)
        w_bomba_real = w_bomba_ideal / eta_bomba
        h2 = h1 + w_bomba_real

        T2 = PropsSI('T', 'P', P_alta, 'H', h2, fluido)
        s2 = PropsSI('S', 'P', P_alta, 'H', h2, fluido)

        # estado 3: entrada da turbina (vapor superaquecido, valor medido)
        h3 = PropsSI('H', 'P', P_alta, 'T', T_3, fluido)
        s3 = PropsSI('S', 'P', P_alta, 'T', T_3, fluido)

        # estado 4: saída da turbina
        # expansão ideal (isentrópica): s4s = s3
        h4_ideal = PropsSI('H', 'P', P_baixa, 'S', s3, fluido)
        T4_ideal = PropsSI('T', 'P', P_baixa, 'S', s3, fluido)
        w_turb_ideal = h3 - h4_ideal
        w_turb_real = w_turb_ideal * eta_turb
        h4 = h3 - w_turb_real

        T4 = PropsSI('T', 'P', P_baixa, 'H', h4, fluido)
        s4 = PropsSI('S', 'P', P_baixa, 'H', h4, fluido)
        # título do vapor: de 0 a 1 se for mistura; CoolProp devolve -1 para fase única
        titulo_4 = PropsSI('Q', 'P', P_baixa, 'H', h4, fluido)
        h_vap_sat = PropsSI('H', 'P', P_baixa, 'Q', 1, fluido)

    except Exception as e:
        print("\nERRO: os dados informados estão fora do domínio de cálculo do CoolProp.")
        print(f"detalhe do erro: {e}")
        return

    # balanço de energia
    w_liquido = (w_turb_real - w_bomba_real) / 1000   # kJ/kg
    q_in = (h3 - h2) / 1000                           # calor fornecido na caldeira, kJ/kg
    q_out = (h4 - h1) / 1000                          # calor rejeitado no condensador, kJ/kg
    eficiencia = (w_liquido / q_in) * 100
    eficiencia_carnot = (1 - T1 / T_3) * 100          # limite: mesmas temperaturas extremas

    print("\n" + "=" * 55)
    print("      estados do ciclo")
    print("=" * 55)
    print(f"{'estado':<8}{'T (°C)':>10}{'P (bar)':>10}{'h (kJ/kg)':>12}{'s (kJ/kg.K)':>14}")
    estados = [
        (1, T1, P_baixa, h1, s1),
        (2, T2, P_alta, h2, s2),
        (3, T_3, P_alta, h3, s3),
        (4, T4, P_baixa, h4, s4),
    ]
    for num, T, P, h, s in estados:
        print(f"{num:<8}{T - 273.15:>10.2f}{P / 1e5:>10.2f}{h / 1000:>12.2f}{s / 1000:>14.4f}")

    print("\n" + "=" * 55)
    print("      desempenho do ciclo")
    print("=" * 55)
    print(f"trabalho específico da bomba:   {w_bomba_real / 1000:8.2f} kJ/kg")
    print(f"trabalho específico da turbina: {w_turb_real / 1000:8.2f} kJ/kg")
    print(f"trabalho líquido específico:    {w_liquido:8.2f} kJ/kg")
    print(f"calor fornecido na caldeira:    {q_in:8.2f} kJ/kg")
    print(f"calor rejeitado no condensador: {q_out:8.2f} kJ/kg")
    print("-" * 55)
    print(f"EFICIÊNCIA TÉRMICA GLOBAL:      {eficiencia:8.2f} %")
    print(f"eficiência de Carnot (limite):  {eficiencia_carnot:8.2f} %")

    if 0 <= titulo_4 <= 1:
        print(f"título de vapor na turbina:     {titulo_4 * 100:8.2f} %")
        if titulo_4 < 0.88:
            print(">>> ALERTA: título baixo (< 88%), risco de erosão nas pás da turbina.")
    elif h4 >= h_vap_sat:
        print("saída da turbina:               vapor superaquecido (sem líquido)")

    # diagrama T-s
    print("\ndiagrama T-s...")

    # curva de saturação (domo da água), do ponto triplo até logo abaixo do ponto crítico
    T_crit = PropsSI('Tcrit', fluido)
    T_domo = np.linspace(273.16, T_crit - 0.01, 200)
    s_liq = [PropsSI('S', 'T', T, 'Q', 0, fluido) / 1000 for T in T_domo]
    s_vap = [PropsSI('S', 'T', T, 'Q', 1, fluido) / 1000 for T in T_domo]

    plt.figure(figsize=(10, 6))
    plt.plot(s_liq, T_domo - 273.15, 'k-', linewidth=1.5, label='Curva de Saturação')
    plt.plot(s_vap, T_domo - 273.15, 'k-', linewidth=1.5)
    # fecha a pequena abertura no topo do domo (ponto crítico)
    plt.plot([s_liq[-1], s_vap[-1]], [T_domo[-1] - 273.15] * 2, 'k-', linewidth=1.5)

    # 1 -> 2: bomba (quase vertical, a entropia varia muito pouco)
    plt.plot([s1 / 1000, s2 / 1000], [T1 - 273.15, T2 - 273.15], 'r-', linewidth=2,
             label='Ciclo Real')

    # 2 -> 3: caldeira (isobárica: aquece o líquido, ferve e superaquece o vapor)
    s_cald, T_cald = linha_isobarica(fluido, P_alta, h2, h3)
    plt.plot(s_cald, T_cald, 'r-', linewidth=2)

    # 3 -> 4: turbina real (irreversível, por isso tracejada)
    plt.plot([s3 / 1000, s4 / 1000], [T_3 - 273.15, T4 - 273.15], 'r--', linewidth=2)

    # 4 -> 1: condensador (isobárica)
    s_cond, T_cond = linha_isobarica(fluido, P_baixa, h4, h1)
    plt.plot(s_cond, T_cond, 'r-', linewidth=2)

    # expansão ideal 3 -> 4s (isentrópica) para comparação
    plt.plot([s3 / 1000, s3 / 1000], [T_3 - 273.15, T4_ideal - 273.15], '--', color='gray',
             linewidth=1.5, label='Expansão ideal (isentrópica)')

    # pontos dos estados e rótulos
    entropias = [s1 / 1000, s2 / 1000, s3 / 1000, s4 / 1000]
    temperaturas = [T1 - 273.15, T2 - 273.15, T_3 - 273.15, T4 - 273.15]
    plt.plot(entropias, temperaturas, 'ro', markersize=8)

    rotulos = ['1 (Entrada da Bomba)', '2 (Saída da Bomba)',
               '3 (Entrada da Turbina)', '4 (Entrada do Condensador)']
    # os estados 1 e 2 quase coincidem (a bomba pouco altera T e s),
    # então os rótulos ficam em lados diferentes para não se sobreporem
    deslocamentos = [(12, -20, 'left'), (-12, 12, 'right'), (12, -10, 'left'), (12, -10, 'left')]
    for i, txt in enumerate(rotulos):
        dx, dy, alinhamento = deslocamentos[i]
        plt.annotate(txt, (entropias[i], temperaturas[i]), textcoords="offset points",
                     xytext=(dx, dy), ha=alinhamento)

    plt.title('Diagrama T-s do Ciclo de Rankine Real (H2O)', fontsize=14, fontweight='bold')
    plt.xlabel('Entropia, s (kJ/kg.K)', fontsize=12)
    plt.ylabel('Temperatura, T (°C)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='upper left')

    # salva o arquivo na pasta e exibe na tela
    plt.savefig("rankine_diagrama_ts.png", dpi=300)
    plt.show()


# ponto de entrada do script
if __name__ == "__main__":
    simular_rankine()
