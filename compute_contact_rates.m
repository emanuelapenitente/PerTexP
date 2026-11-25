function rates = compute_contact_rates()
%COMPUTE_CONTACT_RATES Compute weekly PerTexP contact rates from Mistry data.
%
%   rates = compute_contact_rates()
%
% This function has NO input arguments. It expects the original Mistry
% CSV files to be stored in
%
%   PerTexP/data/mistry_contact_data/
%
% with the filenames
%
%   Italy_country_level_M_overall_contact_matrix_85.csv
%   Italy_country_level_age_distribution_85.csv
%
% The Mistry overall contact matrix contains DAILY contact rates in
% one-year age classes:
%
%   0, 1, 2, ..., 83, 84+
%
% PerTexP uses two age groups:
%
%   group 1: infants, age [0,1)
%   group 2: non-infants, age >= 1 year.
%
% The function:
%   1. reads the original Mistry CSV files;
%   2. converts the 85x85 contact matrix from daily to weekly rates;
%   3. aggregates it into the two PerTexP age groups;
%   4. returns the weekly contact rates C11, C12, C21, C22.
%
% OUTPUT
% ------
% rates.C11_week : contacts of an infant with infants
% rates.C12_week : contacts of an infant with non-infants
% rates.C21_week : contacts of a non-infant with infants
% rates.C22_week : contacts of a non-infant with non-infants
% rates.C_week    : 2x2 weekly contact matrix


    %% Locate repository and data directory
    thisFile = mfilename('fullpath');
    repoDir  = fileparts(thisFile);

    dataDir = fullfile(repoDir, 'data', 'mistry_contact_data');

    contactFile = fullfile(dataDir, ...
        'Italy_country_level_M_overall_contact_matrix_85.csv');

    populationFile = fullfile(dataDir, ...
        'Italy_country_level_age_distribution_85.csv');

    %% Check input files
    if ~isfile(contactFile)
        error(['Mistry contact matrix not found:\n%s\n\n', ...
               'Expected repository layout:\n', ...
               'PerTexP/data/mistry_contact_data/'], contactFile);
    end

    if ~isfile(populationFile)
        error(['Mistry age distribution not found:\n%s\n\n', ...
               'Expected repository layout:\n', ...
               'PerTexP/data/mistry_contact_data/'], populationFile);
    end

    %% Read original Mistry data
    M_day = readmatrix(contactFile);
    age_distribution = readmatrix(populationFile);

    %% Validate contact matrix
    if ~isequal(size(M_day), [85, 85])
        error('Expected an 85x85 contact matrix, but found %dx%d.', ...
              size(M_day,1), size(M_day,2));
    end

    if any(~isfinite(M_day(:))) || any(M_day(:) < 0)
        error('The contact matrix contains invalid values.');
    end

    %% Validate age distribution
    if size(age_distribution,1) ~= 85 || size(age_distribution,2) < 2
        error(['Expected the age-distribution CSV to contain 85 rows ', ...
               'and at least two columns [age, population].']);
    end

    ages = age_distribution(:,1);
    population = age_distribution(:,2);

    if any(abs(ages - (0:84)') > 1e-12)
        error(['The first column of the age-distribution file must ', ...
               'contain age indices 0,1,...,84, where 84 denotes 84+.']);
    end

    if any(~isfinite(population)) || any(population <= 0)
        error('The population column contains invalid values.');
    end

    %% Convert daily contact rates to weekly contact rates
    days_per_week = 7;
    M_week = days_per_week * M_day;

    %% Population sizes in the Mistry age distribution
    %
    % MATLAB index 1 = age 0
    % MATLAB indices 2:85 = ages 1,...,83,84+
    N1_mistry = population(1);
    N2_mistry = sum(population(2:end));

    %% Aggregate DAILY rates into the two PerTexP age groups
    %
    % C11: one age-0 infant with age-0 infants
    C11_day = M_day(1,1);

    % C12: one age-0 infant with all non-infant ages
    C12_day = sum(M_day(1,2:end));

    % C21: average contacts of a non-infant with age-0 infants.
    C21_day = sum(population(2:end) .* M_day(2:end,1)) / N2_mistry;

    % C22: for each non-infant age, first sum contacts with all
    contacts_with_noninfants_day = sum(M_day(2:end,2:end), 2);

    C22_day = sum(population(2:end) .* contacts_with_noninfants_day) ...
              / N2_mistry;

    %% Weekly aggregated rates
    C11_week = days_per_week * C11_day;
    C12_week = days_per_week * C12_day;
    C21_week = days_per_week * C21_day;
    C22_week = days_per_week * C22_day;

    C_day = [C11_day, C12_day; ...
             C21_day, C22_day];

    C_week = [C11_week, C12_week; ...
              C21_week, C22_week];

    %% Reciprocity check for the original Mistry population
    %
    % A reciprocal contact matrix should satisfy
    %
    %   N1*C12 = N2*C21.
    reciprocity_day = N1_mistry*C12_day - N2_mistry*C21_day;
    reciprocity_week = N1_mistry*C12_week - N2_mistry*C21_week;

    tol = 1e-8 * max(1, abs(N1_mistry*C12_week));

    if abs(reciprocity_week) > tol
        warning(['Aggregated contact matrix does not satisfy reciprocity ', ...
                 'within numerical tolerance. Difference = %.6e.'], ...
                 reciprocity_week);
    end

    %% Return results
    rates = struct();

    rates.C11_week = C11_week;
    rates.C12_week = C12_week;
    rates.C21_week = C21_week;
    rates.C22_week = C22_week;
    rates.C_week = C_week;

    rates.C11_day = C11_day;
    rates.C12_day = C12_day;
    rates.C21_day = C21_day;
    rates.C22_day = C22_day;
    rates.C_day = C_day;

    rates.M_day = M_day;
    rates.M_week = M_week;

    rates.N1_mistry = N1_mistry;
    rates.N2_mistry = N2_mistry;

    rates.reciprocity_residual_day = reciprocity_day;
    rates.reciprocity_residual_week = reciprocity_week;

    rates.contact_file = contactFile;
    rates.population_file = populationFile;
    rates.units = 'contacts/person/week';
    rates.source = 'Mistry et al., Italy overall contact matrix, 85 age classes';

    %% Display weekly rates
    fprintf('\n====================================================\n');
    fprintf('Weekly PerTexP contact rates from Mistry et al.\n');
    fprintf('====================================================\n');
    fprintf('C11 = %.12f contacts/person/week\n', C11_week);
    fprintf('C12 = %.12f contacts/person/week\n', C12_week);
    fprintf('C21 = %.12f contacts/person/week\n', C21_week);
    fprintf('C22 = %.12f contacts/person/week\n', C22_week);

    fprintf('\nWeekly 2x2 contact matrix:\n');
    disp(C_week);

    fprintf('Reciprocity residual N1*C12 - N2*C21 = %.6e\n', ...
            reciprocity_week);
end
