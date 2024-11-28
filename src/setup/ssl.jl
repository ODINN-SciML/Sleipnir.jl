# Set parameters of CondaPkg with correct SSL version
using PreferenceTools

SSL_version = "3.4.0"

PreferenceTools.add("CondaPkg", "openssl_version" => SSL_version)

@info("Currently using SSL version $SSL_version of SSL for CondaPkg.")