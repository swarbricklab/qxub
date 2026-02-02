# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 3.x.x   | :white_check_mark: |
| 2.x.x   | :x:                |
| < 2.0   | :x:                |

## Reporting a Vulnerability

We take the security of qxub seriously. If you believe you have found a security vulnerability, please report it responsibly.

### How to Report

**Please do NOT report security vulnerabilities through public GitHub issues.**

Instead, please report them to us privately:

1. **Email**: Send details to the repository maintainers
2. **GitHub Security Advisory**: Use GitHub's private vulnerability reporting feature
3. **Direct Contact**: Contact maintainers directly through institutional channels

### What to Include

Please include the following information in your report:

- **Description** of the vulnerability
- **Steps to reproduce** the issue
- **Potential impact** assessment
- **Suggested fix** (if you have one)
- **Your contact information** for follow-up

### Response Timeline

- **Initial Response**: Within 48 hours
- **Status Update**: Within 1 week
- **Resolution**: Target within 30 days (varies by complexity)

### Disclosure Policy

- We will work with you to understand and resolve the issue
- We will provide credit for your discovery (unless you prefer to remain anonymous)
- We will coordinate public disclosure after the issue is resolved
- We follow responsible disclosure practices

## Security Best Practices

### For Contributors

- **Never commit secrets**: Use environment variables and secure secret management
- **Review dependencies**: Keep dependencies up to date and review security advisories
- **Follow least privilege**: Request minimal permissions in workflows and code
- **Validate inputs**: Sanitize and validate all external inputs
- **Use secure communication**: Always use HTTPS/SSH for remote connections

### For Users

- **Keep qxub updated**: Use the latest version for security fixes
- **Secure your credentials**: Use SSH keys, avoid passwords, rotate credentials regularly
- **Review configurations**: Audit your qxub configurations and remote access settings
- **Monitor job logs**: Check for suspicious activity in job execution logs
- **Use dedicated accounts**: Consider dedicated service accounts for CI/automated systems

## Known Security Considerations

### Remote Execution

- **SSH Key Management**: qxub requires SSH access to remote systems
- **Command Injection**: Commands are executed on remote systems
- **File Transfer**: Files may be transferred between local and remote systems
- **Credential Storage**: SSH keys and configurations require secure storage

### HPC Environment Integration

- **PBS/Slurm Integration**: Job scripts are submitted to schedulers
- **Environment Modules**: System modules and conda environments are activated
- **Shared Storage**: Jobs may access shared file systems
- **Resource Allocation**: Jobs request compute resources from shared pools

## Security Updates

Security updates will be released as patch versions and communicated through:

- **GitHub Releases** with security labels
- **Repository Security Tab** for security advisories
- **CHANGELOG.md** with security fix annotations

## Contact

For security-related questions or concerns, please contact the repository maintainers through appropriate institutional channels.
