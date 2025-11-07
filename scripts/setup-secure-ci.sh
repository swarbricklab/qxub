#!/bin/bash
# Setup script for qxub secure CI with GCP Secret Manager
# Run this script to create the required secrets in GCP

set -e

PROJECT_ID="${PROJECT_ID:-nci-automation}"
SERVICE_ACCOUNT="grunner-ci-sa@${PROJECT_ID}.iam.gserviceaccount.com"

echo "🔐 Setting up qxub secure CI with GCP Secret Manager"
echo "Project ID: $PROJECT_ID"
echo "Service Account: $SERVICE_ACCOUNT"
echo

# Check if gcloud is authenticated
if ! gcloud auth list --filter=status:ACTIVE --format="value(account)" | grep -q .; then
    echo "❌ Please authenticate with gcloud first:"
    echo "   gcloud auth login"
    exit 1
fi

# Check if private key file exists
if [[ ! -f "gadi_private_key" ]]; then
    echo "❌ Please place your Gadi SSH private key in 'gadi_private_key'"
    echo "   This should be the private key that corresponds to your Gadi public key"
    exit 1
fi

# Generate known_hosts for Gadi
echo "🔍 Generating SSH known_hosts for Gadi..."
ssh-keyscan -H gadi.nci.org.au > gadi_known_hosts 2>/dev/null

# Create or update secrets
echo "📝 Creating/updating GCP secrets..."

# Create private key secret
echo "   Creating gadi_ssh_private_key secret..."
if gcloud secrets describe gadi_ssh_private_key >/dev/null 2>&1; then
    echo "   Secret exists, creating new version..."
    gcloud secrets versions add gadi_ssh_private_key --data-file=gadi_private_key
else
    echo "   Creating new secret..."
    gcloud secrets create gadi_ssh_private_key --data-file=gadi_private_key
fi

# Create known_hosts secret
echo "   Creating gadi_ssh_known_hosts secret..."
if gcloud secrets describe gadi_ssh_known_hosts >/dev/null 2>&1; then
    echo "   Secret exists, creating new version..."
    gcloud secrets versions add gadi_ssh_known_hosts --data-file=gadi_known_hosts
else
    echo "   Creating new secret..."
    gcloud secrets create gadi_ssh_known_hosts --data-file=gadi_known_hosts
fi

# Grant access to service account
echo "🔑 Granting access to CI service account..."
for secret in gadi_ssh_private_key gadi_ssh_known_hosts; do
    echo "   Granting access to $secret..."
    gcloud secrets add-iam-policy-binding $secret \
        --member="serviceAccount:$SERVICE_ACCOUNT" \
        --role="roles/secretmanager.secretAccessor" >/dev/null
done

# Cleanup temporary files
rm -f gadi_known_hosts

echo "✅ Setup complete!"
echo
echo "📋 Next steps:"
echo "1. Commit and push the new secure workflow:"
echo "   git add .github/workflows/test-remote-execution-secure.yml docs/secure-ci-setup.md"
echo "   git commit -m 'feat: Add secure CI workflow with GCP Secret Manager'"
echo "   git push"
echo
echo "2. Configure GitHub repository secrets (if not already set):"
echo "   PROJECT_ID: $PROJECT_ID"
echo "   PROJECT_NUMBER: (your GCP project number)"
echo "   NCI_USERNAME: (your Gadi username)"
echo "   NCI_PROJECT: (your NCI project code)"
echo
echo "3. Configure GitHub repository variables (optional):"
echo "   ENABLE_ACTUAL_REMOTE_TEST: 'true' (for real job testing)"
echo
echo "4. Test the workflow:"
echo "   gh workflow run test-remote-execution-secure.yml --repo swarbricklab/qxub"
echo
echo "🔒 Your SSH keys are now securely stored in GCP Secret Manager!"
