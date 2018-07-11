import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib.backends.backend_pdf import PdfPages
from bokeh.palettes import mpl
from bokeh.plotting import figure, show, output_file
from bokeh.layouts import column
from itertools import combinations
import os
import scipy.constants as const

class Deconvolution_1D:
    def __init__(self, mixture, endmember, path, fname):
        self.mixture = mixture
        self.endmember = endmember
        self.path = path
        self.fname = fname

    def planck(self, wave,Temp):
        # Note: in SI units
        pi = const.pi
        c = const.c
        h = const.h
        k = const.k
        # planck's law: Thermal radiation power of a black body per unit area of radiating surface per unit of solid angle and per unit frequency
        numer = 2. * h * (c) ** 2 / wave ** 5
        denom = exp(h * c / (k * temp * wave)) - 1.
        return numer / denom

    def interpolate_lambda(self, mix_wavelength):
        """ This class resizes the endmember wavelength equal to that of the mixture wavelength.
        It carries two arguments (endmember dictionary, wavelength vector of the mixture)"""
        interpolated = {}

        for k in self.endmember.keys():
            xendm = self.endmember[k][:, 0].copy()
            yendm = self.endmember[k][:, 1].copy()
            f = interpolate.interp1d(xendm, yendm)
            interpolated[k] = np.column_stack((mix_wavelength, f(mix_wavelength)))
        return interpolated
    
    def scale(mix,endm):
        high = np.max(mix)
        low = np.min(mix)
        mins = np.min(endm.values, axis=0)
        maxs = np.max(endm.values, axis=0)
        rnge = maxs - mins
        scaled_endm = pd.DataFrame((high - (((high - low) * (maxs - endm)) / rnge)), columns=list(endm))
        return mix, scaled_endm

    def Bokehplot(self, key, wavelength, X_mixture, X_modelled, RMS, chi_en, endm_en):
        # bpkeh plot
        fig = figure(title=key + ' : \n RMS = %f' % RMS, x_axis_label="Wavelength (microns)", y_axis_label="Emissivity", plot_width=1000, plot_height=500)
        fig.title.text_font = "times"
        fig.title.text_font_style = "bold"
        fig.title.text_font_size = '18pt'
        fig.axis.axis_label_text_font_size = '12pt'
        fig.axis.major_label_text_font_size = '12pt'
        fig.axis.major_label_text_font_style = 'bold'
        fig.axis.axis_label_text_font_style = 'bold'
        # fig.xgrid.grid_line_color = None
        # fig.ygrid.grid_line_color = None

        colr_main = ['red', 'green']
        leg_main = ['Telescope Spectra', 'Modelled Spectra']
        for (colr, leg, x, y) in zip(colr_main, leg_main, [wavelength, wavelength], [X_mixture, X_modelled]):
            my_plot = fig.line(x, y, line_width=3, color=colr, legend=leg)

        colr_list = mpl["Viridis"][4][:len(endm_en.keys())]
        cl = 0
        for i in endm_en.keys():
            label = i + ': modal % =' + str(np.array(chi_en[i]))
            my_plot = fig.line(wavelength, np.asarray(endm_en[i]), line_dash=[6, 6], line_width=1,
                               color=colr_list[cl], legend=label)
            cl += 1
        fig.legend.background_fill_alpha = 0
        return fig

    def linear_unmixing(self):
        # adapted from Ramsay and Christensen 1998
        # improved by applying combination of endmembers for each iteration
        # added blackbody

        output_file(os.path.join(self.path, ''.join(self.fname+".html")))

        bkp = [None] * len(self.mixture)
        bk = 0
        for key in self.mixture.keys():
            print(key)
            wavelength = self.mixture[key][:, 0]
            m = len(wavelength)  # (m) wavelength channels
            BB = {"Blackbody": np.ones(m)}
            endm_dict_full = self.interpolate_lambda(wavelength)  # regrid endmembers to actual mixture wavelength
            endm_dict = {j: endm_dict_full[j][:, 1] for j in endm_dict_full.keys()}  # endmember emissivity spectrum

            # endm_iter = [{**{j: endm_dict[j] for j in i},**BB} for i in combinations(endm_dict, 2)]

            endm_iter = [{j: endm_dict[j] for j in i} for i in combinations(endm_dict, 3)]

            ############## initialise first iteration #############
            endm_df = pd.DataFrame(endm_iter[0])  # conversion of dict to Dataframe
            endm_df = endm_df.loc[:, (endm_df > 0).any(axis=0)]  # removing columns (endmembers) that have values <=0

            # ******** scaling to high and low of mixture spectra
            mixture_spectra, endm_df = scale(self.mixture[key][:, 1],endm_df)  # scaling to high and low of mixture spectra
            #
            #
            # # # ******* normalising at last wavelength channel
            # mixture_spectra, endm_df = normalisation.normalise2lambda(self.mixture[key][:, 1], endm_df)

            # mixture_spectra = self.mixture[key][:, 1]
            U = np.matrix(mixture_spectra)  # Unknown mixture spectrum

            while True:  # iteration of chi until chi >0
                X = np.matrix(endm_df.values)  # endmember matrix
                chi = ((((X.T) * X).I) * (X.T)) * (U.T)  # chi equation
                chi = chi / chi.sum()  # column vector summed to unity

                if np.all(chi > 0):
                    chi_df = pd.DataFrame(chi.T, columns=list(endm_df))
                    X = endm_df[list(chi_df)]
                    break

                chi_df = pd.DataFrame(chi.T, columns=list(endm_df))  # chi to dataframe
                chi_df = chi_df.loc[:, (chi_df > 0).any(axis=0)]  # flagging only positive chi
                endm_df = endm_df[list(chi_df)]  # extracting only corresponding positive chi

            X_normalised = pd.DataFrame(np.asarray(chi_df) * np.asarray(X), columns=list(chi_df))

            X_modelled = np.sum(np.asarray(X_normalised), axis=1)
            X_mixture = mixture_spectra  # Oliver -> Unknown mixture spectrum

            endm_en = endm_df.copy()
            chi_en = chi_df.copy()

            delta = np.array(U) - X_modelled
            RMS = np.sqrt(np.sum(delta ** 2) / m)

            #### creating dict to write to file
            chi_dict = chi_df.to_dict('records')[0]
            text_dict = {i: [] for i in endm_dict_full.keys()}  # writing to text file
            for i in text_dict.keys():
                if i in chi_dict.keys():
                    text_dict[i].append(chi_dict[i])
                else:
                    text_dict[i].append(0.0)
            text_dict["RMS"] = [RMS]
            ####


            bkp[bk] = self.Bokehplot(key, wavelength, X_mixture, X_modelled, RMS, chi_en, endm_en)

            #####################################################

            for en in endm_iter[1:]:
                endm_df = pd.DataFrame(en)  # conversion of dict to Dataframe
                endm_df = endm_df.loc[:, (endm_df > 0).any(axis=0)]  # removing columns (endmembers) that have values <=0

                # ******** scaling to high and low of mixture spectra
                mixture_spectra, endm_df = normalisation.scale(self.mixture[key][:, 1], endm_df)  # scaling to high and low of mixture spe


                # ******* normalising at last wavelength channel
                # mixture_spectra, endm_df = normalisation.normalise2lambda(self.mixture[key][:, 1],endm_df)

                # mixture_spectra = self.mixture[key][:, 1]
                U = np.matrix(mixture_spectra)  # Unknown mixture spectrum

                while True:  # iteration of chi until chi >0
                    X = np.matrix(endm_df.values)  # endmember matrix
                    chi = ((((X.T) * X).I) * (X.T)) * (U.T)  # chi equation
                    chi = chi / chi.sum()  # column vector summed to unity

                    if np.all(chi > 0):
                        chi_df = pd.DataFrame(chi.T, columns=list(endm_df))
                        X = endm_df[list(chi_df)]
                        break

                    chi_df = pd.DataFrame(chi.T, columns=list(endm_df))  # chi to dataframe
                    chi_df = chi_df.loc[:, (chi_df > 0).any(axis=0)]  # flagging only positive chi
                    endm_df = endm_df[list(chi_df)]  # extracting only corresponding positive chi

                X_normalised = pd.DataFrame(np.asarray(chi_df) * np.asarray(X), columns=list(chi_df))

                X_modelled_1 = np.sum(np.asarray(X_normalised), axis=1)
                X_mixture_1 = mixture_spectra  # Oliver -> Unknown mixture spectrum

                endm_en_1 = endm_df.copy()
                chi_en_1 = chi_df.copy()

                delta = np.array(U) - X_modelled_1
                RMS_1 = np.sqrt(np.sum(delta ** 2) / m)

                #### dict to write to file
                chi_dict = chi_df.to_dict('records')[0]
                for i in text_dict.keys():
                    if i in chi_dict.keys():
                        text_dict[i].append(chi_dict[i])
                    elif i is "RMS":
                        text_dict[i].append(RMS_1)
                    else:
                        text_dict[i].append(0.0)
                ####

                if RMS < RMS_1:
                    pass

                else:
                    X_modelled = X_modelled_1.copy()
                    X_mixture = X_mixture_1.copy()
                    endm_en = endm_en_1.copy()
                    chi_en = chi_en_1.copy()
                    RMS = RMS_1.copy()
                    bkp[bk] = self.Bokehplot(key, wavelength, X_mixture, X_modelled, RMS, chi_en, endm_en)
            bk += 1
            pd.DataFrame(text_dict).to_csv(os.path.join(self.path, ''.join(self.fname+"_%s.csv") % key), header=True)

        show(column(bkp))