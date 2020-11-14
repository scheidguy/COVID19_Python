
''' Creates plots of USA's COVID cases, because why not
###
### Ideas for future enhancements:
### 1. More sophisticated model of Prevalence Ratio that makes it a
###    function of both test positivity rate and number of tests per capita
###    and allows the [a b c] parameters to vary over time:
###    https://covid19-projections.com/estimating-true-infections/
###
### 2. County level and Nationwide plots that incorporate Prevalance Ratio
###
### 3. Incorporating estimates of unreported COVID deaths:
###    https://weinbergerlab.github.io/excess_pi_covid/
###
### 4. Quantifying the amount of correlation between shutdowns/mask
###    mandates/etc and number of new infections and/or R_t. This seems
###    hard...and FUN!
'''


import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import binom
import matplotlib.pyplot as plt 


makeFigsFlag = False
makeCountyFigsFlag = False
value = ''
while value != 'y' and value != 'n':
    value = input("Make County Figs?\n")
if value == 'y':
    makeCountyFigsFlag = True
    makeFigsFlag = True
if not makeCountyFigsFlag:
    value = ''
    while value != 'y' and value != 'n':
        value = input("Make State Figs?\n")
    if value == 'y':
        makeFigsFlag = True


# Load some data. Can find new data from these websites:
url = 'https://usafactsstatic.blob.core.windows.net/public/data/covid-19/'
url2 = 'https://api.covidtracking.com/v1/states/'
url3 = 'https://d14wlfuexuxgcm.cloudfront.net/covid/'
popFile = 'covid_county_population_usafacts.csv'
caseFile = 'covid_confirmed_usafacts.csv'
deathFile = 'covid_deaths_usafacts.csv'
testFile = 'daily.csv'
RTfile = 'rt.csv'
popdf = pd.read_csv(url + popFile)
population = popdf.iloc[:, 3]
counties = popdf.iloc[:, 1]
states = popdf.iloc[:, 2]
casedf = pd.read_csv(url + caseFile)
cases = casedf.iloc[:, 4:]
deathdf = pd.read_csv(url + deathFile)
deaths = deathdf.iloc[:, 4:]
testdf = pd.read_csv(url2 + testFile)
testsStates = testdf.loc[:, 'state']
testsStates = np.array(testsStates[::-1])  # Want oldest records first
tests = testdf[['positive', 'negative', 'date']]  # number of pos and neg tests
tests = tests[::-1]  # oldest records first
RTdf = pd.read_csv(url3 + RTfile)
RTliveStates = RTdf['region']
RTliveValues = RTdf['mean']
RTliveValues20 = RTdf['lower_80']
RTliveValues80 = RTdf['upper_80']


infectionLength = 15  # Avg infection length suggested by covid-projections
daysPerMonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]  # Leap Day!
offsetDays = 21  # Because first date in some spreadsheets is January 22nd


# Have to address the weird New Jersey artifact on July 17th of excess
# deaths. Choosing to just ignore it for now and copy the data from July
# 16th instead.
theDiff = deaths.iloc[1807:1829, 157] - deaths.iloc[1807:1829, 156]
row = 0
for diff in theDiff:
    deaths.iloc[1807+row, 157:] = deaths.iloc[1807+row, 157:] - diff
    row += 1


# Calculate per million multiplier and fix the divide by zeros
multiplier = 10**6 / population
multiplier = multiplier.replace([np.inf, -np.inf], 0)


# Initialize arrays to store all the data and a couple counters
dimensions = cases.shape
casesPerDay = np.zeros(dimensions)
deathsPerDay = np.zeros(dimensions)
cases7day = np.zeros(dimensions)
deaths7day = np.zeros(dimensions)
cases7dayPerMil = np.zeros(dimensions)
deaths7dayPerMil = np.zeros(dimensions)
casesPerDayStatewide = np.zeros((51, dimensions[1]))
deathsPerDayStatewide = np.zeros((51, dimensions[1]))
cases7dayStatewide = np.zeros((51, dimensions[1]))
deaths7dayStatewide = np.zeros((51, dimensions[1]))
cases7dayPerMilStatewide = np.zeros((51, dimensions[1]))
deaths7dayPerMilStatewide = np.zeros((51, dimensions[1]))
deathsStatewideCumul = np.zeros((51, dimensions[1]))
estimatedPercentInfected = np.zeros((51, 1))
previousState = states[0]
stateOrder = []
stateBeginIndex = 0
currState = -1


# Create 7 day averages for cases and deaths per million population
# Have to loop through each county
for currCounty in range(dimensions[0]):
    # cases/deaths are both cumulative, so need to figure out per day stats
    casesPerDay[currCounty, :] = cases.iloc[currCounty, :].diff()
    casesPerDay[currCounty, 0] = cases.iloc[currCounty, 0]
    deathsPerDay[currCounty, :] = deaths.iloc[currCounty, :].diff()
    deathsPerDay[currCounty, 0] = deaths.iloc[currCounty, 0]

    # Have to loop through each day
    for currDay in range(dimensions[1]):
        # Average everything if haven't hit 7th day yet
        if currDay < 6:
            cases7day[currCounty, currDay] = sum(casesPerDay[currCounty, 0:currDay+1] / (currDay+1))
            deaths7day[currCounty, currDay] = sum(deathsPerDay[currCounty, 0:currDay+1] / (currDay+1))        
        else:   # otherwise just average most recent 7 days
            cases7day[currCounty, currDay] = sum(casesPerDay[currCounty, currDay-6:currDay+1] / 7)
            deaths7day[currCounty, currDay] = sum(deathsPerDay[currCounty, currDay-6:currDay+1] / 7)


    # Create statewide plot if we've iterated to a new state
    newState = False
    if previousState != states[currCounty]:
        newState = True
    if newState or currCounty == dimensions[0]-1:
        currState += 1
        stateOrder.append(previousState)
        if stateBeginIndex != (currCounty-1):  # This means we aren't DC
            casesPerDayStatewide[currState, :] = sum(casesPerDay[stateBeginIndex:currCounty, :])
            deathsPerDayStatewide[currState, :] = sum(deathsPerDay[stateBeginIndex:currCounty, :])
            cases7dayStatewide[currState, :] = sum(cases7day[stateBeginIndex:currCounty, :])
            deaths7dayStatewide[currState, :] = sum(deaths7day[stateBeginIndex:currCounty, :])
            cases7dayPerMilStatewide[currState, :] = sum(cases7dayPerMil[stateBeginIndex:currCounty, :])
            deaths7dayPerMilStatewide[currState, :] = sum(deaths7dayPerMil[stateBeginIndex:currCounty, :])
            deathsStatewideCumul[currState,:] = deaths.iloc[stateBeginIndex:currCounty, :].sum()
            statePop = sum(population[stateBeginIndex:currCounty])
            multiplierStatewide = 10**6 / statePop
        else:  # Otherwise we're in Washington DC
            casesPerDayStatewide[currState, :] = casesPerDay[stateBeginIndex:currCounty, :]
            deathsPerDayStatewide[currState, :] = deathsPerDay[stateBeginIndex:currCounty, :]
            cases7dayStatewide[currState, :] = cases7day[stateBeginIndex:currCounty, :]
            deaths7dayStatewide[currState, :] = deaths7day[stateBeginIndex:currCounty, :]
            cases7dayPerMilStatewide[currState, :] = cases7dayPerMil[stateBeginIndex:currCounty, :]
            deaths7dayPerMilStatewide[currState, :]= deaths7dayPerMil[stateBeginIndex:currCounty, :]
            deathsStatewideCumul[currState,:] = deaths.iloc[stateBeginIndex:currCounty, :].sum()
            statePop = sum(population[stateBeginIndex:currCounty])
            multiplierStatewide = 10**6 / statePop
    
    
        # Calculate the cumulative cases and deaths per million population
        casesPerMilStatewide = round(multiplierStatewide * sum(casesPerDayStatewide[currState,:]), -2)
        deathPerMilStatewide = round(multiplierStatewide * sum(deathsPerDayStatewide[currState,:]), 1)
    
    
        # Calculate test positivity rate
        testInds = np.where(np.isin(testsStates, previousState))
        currTests = tests.iloc[testInds[0], :]
        # Remove rows without a tally for negative tests
        currTests = currTests.replace(np.nan, 0)
        positivityRate = currTests.iloc[:,0] / (currTests.iloc[:,0] + currTests.iloc[:,1])
        positivityRate = positivityRate.replace(np.nan, 0)  # Replace divide-by-0s
    
        # Have to loop through each day to calculate 7 day average
        positivityRate7day = np.zeros((len(positivityRate), 1))
        dayNumber = np.zeros((len(positivityRate), 1))
        confirmedCases = np.zeros((len(positivityRate), 1))
        for currentDay in range(len(positivityRate)):
            # Average everything if haven't hit 7th day yet
            if currentDay < 6:
                positivityRate7day[currentDay] = sum(positivityRate[0:currentDay+1] / (currentDay+1))
            else:   # otherwise just average most recent 7 days
                positivityRate7day[currentDay] = sum(positivityRate[currentDay-6:currentDay+1] / 7)

        
            # Figure out how many days from January 21st we are:
            monthDay = currTests.iloc[currentDay,2] - 20200000
            numMonths = int(round(monthDay/100))
            if monthDay > 1231:
                sys.exit('Need to adjust code to accept non-2020 years!')
            if numMonths > 1:  # Need to add previous months days
                dayNumber[currentDay] = (monthDay % 100) + sum(daysPerMonth[0:numMonths-1])
            else:
                dayNumber[currentDay] = (monthDay % 100)

        
            # Grab confirmed cases 7 day average for this day and state
            theDay = int(dayNumber[currentDay] - offsetDays)
            if theDay <= np.shape(cases7dayStatewide)[1]:
                confirmedCases[currentDay] = cases7dayStatewide[currState, theDay-1]
            else:  # Don't have new versions of the other spreadsheets
                # Just assume same as previous day
                confirmedCases[currentDay] = cases7dayStatewide[currState, theDay-2]

    
        # Calculate estimate of TRUE number of cases as a function of 7 day
        # average of test positivity rate, as described here:
        # https://covid19-projections.com/estimating-true-infections
        MayInd = np.where(dayNumber > sum(daysPerMonth[0:4]))[0][0]
        a1 = 16
        b1 = 0.5
        c1 = 2.5
        prevalenceRatioEarly = a1 * positivityRate7day[0:MayInd] ** b1 + c1
        a2 = 10
        b2 = 0.4
        c2 = 2.5
        prevalenceRatioLater = a2 * positivityRate7day[MayInd:] ** b2 + c2
        prevalenceRatio = np.concatenate([prevalenceRatioEarly, prevalenceRatioLater])
        estimatedCases = prevalenceRatio * confirmedCases
        estimatedPercentInfected[currState] = round(float(100 * (sum(estimatedCases) / statePop)), 1)
        
        
        # Calculate estimate of percent currently infected
        estimatedCurrentlyInfected = np.zeros((len(estimatedCases), 1))
        for dayIter in range(len(estimatedCases)):
            if dayIter < infectionLength:
                estimatedCurrentlyInfected[dayIter] = 100 * (sum(estimatedCases[0:dayIter+1] / statePop))
            else:
                estimatedCurrentlyInfected[dayIter] = 100 * (sum(estimatedCases[dayIter-infectionLength:dayIter+1] / statePop))
        print(f'{previousState}: {round(float(estimatedCurrentlyInfected[-1]), 2)}% Currently Infected')


        # Grab the R_t data for this state
        casesInds = np.where(np.isin(RTliveStates, previousState))
        RTcases = RTliveValues.iloc[casesInds[0]]
        RTcases20 = RTliveValues20.iloc[casesInds[0]]
        RTcases80 = RTliveValues80.iloc[casesInds[0]]
        RTcasesDays = np.arange(int(dayNumber[-1]-len(RTcases)+1), int(dayNumber[-1]+1))
        RTcasesNow = round(float(RTcases.iloc[-1]), 2)
        # Grab data from May 2020 onwards (once testing ramped up) 
        RTMayInd = np.where(RTcasesDays > sum(daysPerMonth[0:4]))[0][0]
        laterRT = np.array(RTcases[RTMayInd:])
        laterEstimated = estimatedCurrentlyInfected[MayInd:]
        laterPositivity = positivityRate7day[MayInd:]
        laterConfirmedCases = casesPerDayStatewide[currState, 121-offsetDays:]
        # Make sure data ends on same day since drawing from multiple sources
        if len(laterRT) > len(laterEstimated):
            laterRT = np.delete(laterRT, -1)
        elif len(laterEstimated) > len(laterRT):
            laterEstimated = np.delete(laterEstimated, -1)
            laterPositivity = np.delete(laterPositivity, -1)
            laterConfirmedCases = np.delete(laterConfirmedCases, -1)
        # Remove most recent week since that is what we are using as "truth"
        for noUse in range(7):
            laterRT = np.delete(laterRT, -1)
            laterEstimated = np.delete(laterEstimated, -1)
            laterPositivity = np.delete(laterPositivity, -1)
        # Calculate # of confirmed cases in following week
        numCasesNextWeek = np.array([])
        for day in range(len(laterRT)):
            numCasesNextWeek = np.append(numCasesNextWeek, sum(laterConfirmedCases[day+1:day+8]))
        state = currState * np.ones((len(numCasesNextWeek), 1))
        pop = statePop * np.ones((len(numCasesNextWeek), 1))
        days = np.arange(122, 122+len(numCasesNextWeek))
        if currState == 0:
            features = np.column_stack((state, pop, days, numCasesNextWeek, laterRT, laterPositivity, laterEstimated))
        else:
            stateFeats = np.column_stack((state, pop, days, numCasesNextWeek, laterRT, laterPositivity, laterEstimated))
            features = np.concatenate((features, stateFeats))


        if makeFigsFlag:
            # Create the figure with 4 subplots
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(16,9))
            fig.suptitle('COVID-19 Stats for the Great State of ' + previousState, size=30)
    
    
            ax1.plot(np.arange(currDay+1)+offsetDays, casesPerDayStatewide[currState,:], color='black', lineWidth=1, label=f'New Confirmed Cases: {int(casesPerMilStatewide)} per million cumulatively')  # 'Color', [0.9290 0.6940 0.1250]);
            ax1.plot(np.arange(currDay+1)+offsetDays, cases7dayStatewide[currState,:], color='cyan', lineWidth=4, label='7 Day Case Average')  # 'LineWidth', 4, 'Color', [0.9290 0.6940 0.1250]);
            ax1.plot(dayNumber, estimatedCases, lineWidth=4, color='blue', label=f'Estimated 7 Day Avg: {float(estimatedPercentInfected[currState])}% of population have had COVID-19')  # 'LineWidth', 4, 'Color', [0.4940 0.1840 0.5560]);
            ax1.set_xlabel('Days Since 2020 Began', fontsize=14)
            ax1.set_ylabel('New COVID-19 Cases', fontsize=14)
            legend = ax1.legend(loc='upper left', fontsize='small')
    
    
            infectedNow = float(estimatedCurrentlyInfected[-1])
            ax2.set_xlabel('Days Since 2020 Began', fontsize=14)
            ax2.set_ylabel('Estimated % Infected at the Time', color='blue', fontsize=14)
            ax2.plot(dayNumber, estimatedCurrentlyInfected, lineWidth=4, color='blue', label=f'% of Population Infected: Currently {round(infectedNow,2)}%')
            ax2.tick_params(axis='y', labelcolor='blue')
            Ax2 = ax2.twinx()
            Ax2.set_ylabel('# of Deaths', color='red', fontsize=14)
            Ax2.plot(np.arange(currDay+1)+offsetDays, deathsStatewideCumul[currState,:], lineWidth=4, color='red', label=f'Cumulative Deaths: {deathPerMilStatewide} per million so far')
            Ax2.tick_params(axis='y', labelcolor='red')
            ax2.legend(loc=(0.01,0.9), fontsize='small')
            Ax2.legend(loc=(0.01,0.8), fontsize='small')
    
    
            ax3.set_xlabel('Days Since 2020 Began', fontsize=14)
            ax3.set_ylabel('Effective Reproduction Number (R_t)', fontsize=14)
            ax3.plot(RTcasesDays, RTcases, lineWidth=4, color='cyan', label=f'R_t Estimate from Cases with 60% Confidence Range: Currently {RTcasesNow}')
            ax3.fill_between(RTcasesDays, RTcases20, RTcases80, color='cyan', alpha=0.2)
            legend = ax3.legend(loc='upper left', fontsize='small')
        
    
            theCDF = 100 * (1 - binom.cdf(0, range(1,1001), infectedNow/100))
            firstOver90 = np.where(theCDF > 90)[0]
            ax4.set_xlabel(f'# of People Interacted With (Assumes {round(infectedNow,2)} % Infected)')
            ax4.set_ylabel('% Chance You Interacted With Infected Person')
            if len(firstOver90) != 0:  # Want to plot up to 90% chance
                ax4.plot(range(1,firstOver90[0]+2), theCDF[0:firstOver90[0]+1], lineWidth=4)
            else:
                ax4.plot(range(1,1001), theCDF, lineWidth=4)
            # Save jpgs
            fig.savefig('jpgs\\'+ previousState + '\\(STATEWIDE).jpg')
            plt.close(fig)

            
        # Update beginning index for the new state
        stateBeginIndex = currCounty
        previousState = states[currCounty]


    # Calculate the cumulative cases and deaths per million population
    casesPerMil = round(float(multiplier[currCounty] * sum(casesPerDay[currCounty,:])), -2)
    deathPerMil = round(float(multiplier[currCounty] * sum(deathsPerDay[currCounty,:])), 1)
    # Generate and save county plots. Don't create plots if "Unallocated"
    if 'Unallocated' in counties[currCounty] or not makeFigsFlag or not makeCountyFigsFlag:
        continue
    fig = plt.figure()
    plt.plot(np.arange(currDay+1)+offsetDays, casesPerDay[currCounty,:], lineWidth=1, color='black', label=f'New Confirmed Cases: {int(casesPerMil)} per million')
    plt.plot(np.arange(currDay+1)+offsetDays, cases7day[currCounty,:], lineWidth=4, color='cyan', label='7 Day Case Average')
    plt.plot(np.arange(currDay+1)+offsetDays, deaths.iloc[currCounty,:], lineWidth=4, color='red', label=f'Cumulative Deaths: {deathPerMil} per million')
    plt.xlabel('Days Since 2020 Began')
    plt.ylabel('New COVID-19 Cases')
    plt.title(f'COVID-19 Stats For {counties[currCounty]}, {states[currCounty]}')
    plt.legend(loc='upper left', fontsize='small')
    saveString = counties[currCounty].replace(' ', '')
    saveString = saveString.replace('.', '')
    # Don't need to save a County fig for Washington DC
    if states[currCounty] != 'DC':
        fig.savefig('jpgs\\' + states[currCounty] + '\\' + saveString + '.jpg')
    plt.close(fig)


# Print some stats to the command line for the worst impacted states
indices = np.where(estimatedPercentInfected >= 20)[0]
for ii in range(len(indices)):
    print(stateOrder[indices[ii]] + ': ' + str(float(estimatedPercentInfected[indices[ii]])) + '% Total Population Infected')


if makeFigsFlag:  # Create national plots
    multiplierNation = 10**6 / sum(population)
    casesPerMilNation = round(multiplierNation * sum(cases.iloc[:,-1]), -2)
    deathPerMilNation = round(multiplierNation * sum(deaths.iloc[:,-1]))
    fig = plt.figure()
    plt.plot(np.arange(currDay+1)+offsetDays, np.sum(casesPerDay,axis=0), lineWidth=1, color='black', label=f'New Confirmed Cases: {casesPerMilNation} per million')
    plt.plot(np.arange(currDay+1)+offsetDays, deaths.sum(axis=0), lineWidth=4, color='red', label=f'Cumulative Deaths: {deathPerMilNation} per million')
    plt.plot(np.arange(currDay+1)+offsetDays, np.sum(cases7day,axis=0), lineWidth=4, color='cyan', label='7 Day Case Average')
    plt.xlabel('Days Since 2020 Began')
    plt.ylabel('New COVID-19 Cases')
    plt.title('COVID-19 Stats for the United States of America')
    plt.legend(loc='upper left', fontsize='small')
    fig.savefig('jpgs\\NATIONWIDE_CASES.jpg')
    plt.close(fig)


    fig = plt.figure()
    plt.plot(np.arange(currDay+1)+offsetDays, np.sum(deathsPerDay,axis=0), lineWidth=1, color='magenta', label=f'COVID-19 Deaths: {deathPerMilNation} per million')
    plt.plot(np.arange(currDay+1)+offsetDays, np.sum(deaths7day,axis=0), lineWidth=4, color='red', label='7 Day Average of Deaths')
    plt.xlabel('Days Since 2020 Began')
    plt.ylabel('New COVID-19 Deaths')
    plt.title('COVID-19 Stats for the United States of America')
    plt.legend(loc='upper left', fontsize='small')
    fig.savefig('jpgs\\NATIONWIDE_DEATHS.jpg')
    plt.close(fig)
