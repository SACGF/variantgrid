const VcSettings = {};
VcSettings.LAB = 'Lab';
VcSettings.INSTITUTION = 'Institution';
VcSettings.LOGGED_IN_USERS = 'Logged in Users';
VcSettings.PUBLIC = '3rd Party Databases';
VcSettings.LOGGED_IN_USERS_MESSAGE = 'This record is shared to logged in users.';

// variantgrid.com enables 'public' ("3rd Party Databases") so curators can share outward
// to ClinVar / MatchMaker Exchange. (Default deployment omits it - opt-in per deployment.)
VcSettings.SHARE_LEVELS = ['lab', 'organisation', 'logged_in_users', 'public'];
