# Domain Pitfalls: GCP Spot VM Deployment

**Domain:** Docker application deployment to GCP spot/preemptible VMs
**Project:** qm-nmr-calc v2.6 - Google Cloud Spot Deployment
**Researched:** 2026-02-04
**Focus:** Long-running NMR calculations (hours) on preemptible infrastructure

---

## Executive Summary

Deploying Docker-based long-running computations to GCP Spot VMs presents unique challenges around preemption handling, data persistence, and deployment reliability. The critical constraint: **GCP provides only 30 seconds warning before termination** - far less than the hours an NMR calculation takes. This document catalogs pitfalls specific to this deployment pattern with actionable prevention strategies.

**Confidence Level:** HIGH - Based on official GCP documentation and community best practices.

---

## Critical Pitfalls

Mistakes that cause data loss, failed deployments, or significant rework.

### Pitfall 1: Ignoring the 30-Second Preemption Window

**What goes wrong:**
NWChem calculations run for hours but GCP terminates Spot VMs with only 30 seconds notice. Without handling this, in-progress calculations are lost entirely, wasting hours of compute time and potentially leaving corrupted job state.

**Why it happens:**
Developers assume preemption is rare or that spot pricing implies some stability guarantee. GCP's aggressive preemption can occur at any time based on demand. Unlike AWS (2 minutes) or Azure (30 seconds with some flexibility), GCP's 30-second window is hard-enforced.

**Consequences:**
- Multi-hour NMR calculations lost entirely
- Job marked as "running" but worker is gone
- User confusion about job status
- Wasted cloud spend on incomplete work

**Prevention:**
1. **Accept calculation loss for v2.6** - For manual lifecycle management, the simplest approach is to accept that preemption means re-running. Document this limitation clearly.
2. **Implement job state persistence** - Save job progress to persistent disk before each major step (geometry opt, NMR calculation per conformer)
3. **Add shutdown script** - Register a metadata shutdown script that marks jobs as "preempted" rather than "running"
4. **Consider NWChem restart files** - NWChem supports restart via `.db` and `.movecs` files. Configure `permanent_dir` to save these to persistent storage.

**Detection (warning signs):**
- Job shows "running" but no progress updates
- Worker container not responding to health checks
- VM shows TERMINATED status in GCP console

**Phase to address:** Infrastructure setup phase - implement shutdown script and job state handling early.

**Sources:**
- [GCP Spot VMs Documentation](https://cloud.google.com/compute/docs/instances/spot)
- [GCP Shutdown Scripts](https://cloud.google.com/compute/docs/shutdownscript)
- [NWChem Restart Documentation](https://nwchemgit.github.io/Start_Restart.html)

---

### Pitfall 2: Using Local SSD Instead of Persistent Disk for Job Data

**What goes wrong:**
Local SSDs provide better performance but **all data is lost when a Spot VM is preempted**. This includes job inputs, intermediate results, and completed outputs.

**Why it happens:**
Local SSDs are faster and cheaper per GB. The documentation trade-off is easy to miss: "Local SSD storage is not automatically replicated and all data can be lost."

**Consequences:**
- Completed job results vanish on preemption
- Cannot recover partial work
- Users lose access to their results without warning

**Prevention:**
1. **Use Persistent Disk for all job data** - Attach a separate persistent disk mounted at `/data` or similar
2. **Configure Docker volumes** - Map job directories to the persistent disk mount point
3. **Set `--no-auto-delete` for data disk** - Ensure the data disk survives even if VM is deleted

**Detection:**
- After VM restart, `/data` directory is empty
- Docker volumes show no data
- Job files referenced in database don't exist on disk

**Phase to address:** Infrastructure setup phase - disk configuration is foundational.

**Sources:**
- [GCP Persistent Disk Documentation](https://cloud.google.com/compute/docs/disks/persistent-disks)
- [Spot VMs and Storage](https://cloud.google.com/compute/docs/instances/spot)

---

### Pitfall 3: Service Account Missing Artifact Registry Permissions

**What goes wrong:**
VM starts but Docker cannot pull images from Google Artifact Registry (or GHCR). Container fails to start silently or with cryptic authentication errors.

**Why it happens:**
The default Compute Engine service account may not have `roles/artifactregistry.reader`. For GHCR, the VM has no authentication configured at all.

**Consequences:**
- Startup script completes but no containers running
- Error logs show "unauthorized" or "manifest not found"
- VM appears healthy but application is down

**Prevention:**
1. **For GCR/Artifact Registry** - Grant `roles/artifactregistry.reader` to the Compute Engine service account
2. **Configure storage access scope** - Use `storage-ro` or `storage-rw` scope when creating VM
3. **For GHCR public images** - No auth needed for public repos (qm-nmr-calc uses public GHCR)
4. **For GHCR private images** - Configure Docker credential helper or use `docker login`

**Detection:**
- `docker ps` shows no containers
- `docker logs` shows auth errors
- Startup script logs show pull failures

**Phase to address:** Deployment script phase - test image pull before full deployment.

**Sources:**
- [Artifact Registry Authentication](https://cloud.google.com/artifact-registry/docs/docker/authentication)
- [Deploying to Compute Engine from Artifact Registry](https://cloud.google.com/artifact-registry/docs/integrate-compute)

---

### Pitfall 4: Firewall Rules Not Created for Custom VPC

**What goes wrong:**
VM is created in a custom VPC (not "default") but SSH and HTTP/HTTPS access is blocked. Cannot connect to the instance or access the web application.

**Why it happens:**
The default VPC has pre-configured firewall rules allowing SSH (port 22), HTTP (80), and HTTPS (443). Custom VPCs start completely locked down - no default rules exist.

**Consequences:**
- Cannot SSH to debug issues
- Web application unreachable
- Health checks fail (if using load balancer)

**Prevention:**
1. **Use default VPC for simplicity** - Unless there's a specific reason for custom networking
2. **Create explicit firewall rules** - SSH from IAP (35.235.240.0/20), HTTPS from 0.0.0.0/0
3. **Verify firewall rules before deployment** - `gcloud compute firewall-rules list`
4. **Use IAP for SSH** - More secure than exposing port 22 to the internet

**Detection:**
- SSH connection timeout
- `curl` to VM external IP fails
- GCP Console shows no ingress firewall rules

**Phase to address:** Infrastructure setup phase - firewall rules are prerequisite for any access.

**Sources:**
- [GCP Firewall Rules](https://cloud.google.com/firewall/docs/firewalls)
- [SSH Best Practices](https://cloud.google.com/compute/docs/connect/ssh-best-practices/network-access)
- [Troubleshooting SSH](https://cloud.google.com/compute/docs/troubleshooting/troubleshooting-ssh-errors)

---

### Pitfall 5: Docker Containers Not Starting After VM Restart

**What goes wrong:**
After a preemption or manual restart, the VM comes back but Docker containers are not running. Application is unreachable.

**Why it happens:**
Multiple causes:
- Container restart policy not set
- Docker Compose not configured to start on boot
- Startup script relies on deprecated container deployment feature
- Container-Optimized OS cloud-init quirks

**Consequences:**
- Application downtime until manual intervention
- Job queue not processing
- User-facing service unavailable

**Prevention:**
1. **Use `restart: unless-stopped` in docker-compose.yml** - Containers restart automatically
2. **Enable Docker service on boot** - `systemctl enable docker`
3. **Use startup script to run `docker compose up -d`** - Idempotent, works on fresh start and restart
4. **Avoid deprecated container deployment feature** - Use explicit docker commands instead

**Detection:**
- `docker ps` shows no running containers after restart
- `docker ps -a` shows containers in "Exited" state
- Application health endpoint not responding

**Phase to address:** Deployment script phase - ensure startup script handles both fresh and restart cases.

**Sources:**
- [Container Deployment on GCE](https://cloud.google.com/compute/docs/containers/deploying-containers)
- [Docker Restart Policies](https://docs.docker.com/config/containers/start-containers-automatically/)

---

## Moderate Pitfalls

Mistakes that cause delays, confusion, or technical debt.

### Pitfall 6: Using `--preemptible` Instead of `--provisioning-model=SPOT`

**What goes wrong:**
Using the legacy `--preemptible` flag instead of the modern `--provisioning-model=SPOT` flag. Works but loses features and follows deprecated path.

**Why it happens:**
Older documentation and blog posts reference `--preemptible`. Google recommends Spot VMs but maintains backward compatibility.

**Consequences:**
- Limited to 24-hour maximum runtime (preemptible VMs are force-terminated at 24 hours)
- Missing newer features
- Following deprecated pattern

**Prevention:**
```bash
# Correct (Spot VM - no 24-hour limit)
gcloud compute instances create my-vm --provisioning-model=SPOT

# Deprecated (Preemptible - 24-hour limit)
gcloud compute instances create my-vm --preemptible
```

**Detection:**
- VM terminates after exactly 24 hours
- Documentation/scripts reference `--preemptible`

**Phase to address:** Deployment script phase - use correct flags from the start.

**Sources:**
- [Create and Use Spot VMs](https://cloud.google.com/compute/docs/instances/create-use-spot)

---

### Pitfall 7: Not Setting Termination Action to STOP

**What goes wrong:**
Default termination action is STOP, but if accidentally set to DELETE, the VM and its metadata are destroyed on preemption. Cannot restart the same VM.

**Why it happens:**
Misunderstanding the `--instance-termination-action` flag or copying examples that use DELETE for stateless workloads.

**Consequences:**
- VM must be recreated from scratch
- Cannot preserve instance metadata
- Persistent disk association may be lost (if auto-delete enabled)

**Prevention:**
```bash
gcloud compute instances create my-vm \
  --provisioning-model=SPOT \
  --instance-termination-action=STOP
```

**Detection:**
- VM disappears from instance list after preemption
- Cannot find VM to restart

**Phase to address:** Deployment script phase - explicit STOP action.

**Sources:**
- [Spot VM Termination Action](https://cloud.google.com/compute/docs/instances/create-use-spot)

---

### Pitfall 8: Spot VM Quota Consuming Standard Quota

**What goes wrong:**
Spot VMs consume standard CPU/GPU quota if no preemptible quota is requested. May hit quota limits meant for production standard VMs.

**Why it happens:**
Preemptible/Spot quota is separate but optional. If not requested, Spot VMs draw from standard quota.

**Consequences:**
- Cannot create standard VMs when needed
- Unexpected quota errors
- Confusion about resource limits

**Prevention:**
1. **Request preemptible quota** - Separate from standard quota
2. **Monitor quota usage** - `gcloud compute regions describe REGION`
3. **For small deployments** - May not matter, standard quota often sufficient

**Detection:**
- Standard VM creation fails with quota error
- Quota dashboard shows high usage from Spot VMs

**Phase to address:** Infrastructure setup phase - consider quota if running multiple VMs.

**Sources:**
- [Spot VMs and Quota](https://cloud.google.com/compute/docs/instances/spot)

---

### Pitfall 9: Premium OS Image Billing Surprise

**What goes wrong:**
Using a premium OS (RHEL, Windows, SLES) with Spot VMs. The Spot discount applies to compute but NOT to OS licensing.

**Why it happens:**
Assumption that Spot pricing applies to everything. OS licensing is billed separately and has minimum usage increments.

**Consequences:**
- Higher than expected bills
- Minimum usage charges even for brief runs
- No cost benefit from preemption under 1 minute

**Prevention:**
1. **Use Debian, Ubuntu, or Container-Optimized OS** - No premium licensing
2. **qm-nmr-calc recommendation: COS or Ubuntu** - Both work well for Docker workloads
3. **Calculate total cost including OS** - If premium OS required

**Detection:**
- Bill higher than Spot pricing calculator suggested
- Line items for OS licensing

**Phase to address:** Infrastructure setup phase - choose OS image carefully.

**Sources:**
- [Preemptible VM Billing](https://cloud.google.com/compute/docs/instances/preemptible)
- [VM Pricing](https://cloud.google.com/compute/vm-instance-pricing)

---

### Pitfall 10: No Health Checks or Status Monitoring

**What goes wrong:**
VM preempted but no alerting. Users discover downtime hours later when checking job status.

**Why it happens:**
Focus on deployment without considering observability. "It's just a spot VM, I'll check manually."

**Consequences:**
- Extended downtime without awareness
- Delayed job restarts
- Poor user experience

**Prevention:**
1. **Use GCP Operations Suite (Stackdriver)** - Uptime checks, alerting
2. **Simple approach: Check VM status programmatically** - `gcloud compute instances describe`
3. **Application-level health endpoint** - `/health` returns job queue status
4. **Email/Slack notification on preemption** - Via shutdown script

**Detection:**
- Find out about downtime from users
- Jobs stalled for hours without notice

**Phase to address:** Operations/monitoring phase - add after core deployment works.

**Sources:**
- [GCP Monitoring](https://cloud.google.com/monitoring)

---

## Minor Pitfalls

Mistakes that cause annoyance but are easily fixable.

### Pitfall 11: Startup Script Not Logging Output

**What goes wrong:**
Startup script fails but no logs to debug. Blind troubleshooting.

**Why it happens:**
Scripts run but output goes nowhere useful. Default logging may not capture script errors.

**Prevention:**
```bash
#!/bin/bash
exec > >(tee -a /var/log/startup-script.log) 2>&1
echo "Starting deployment at $(date)"
# ... rest of script
```

**Detection:**
- Cannot find startup script output
- `/var/log/messages` or journal has no useful info

**Phase to address:** Deployment script phase - add logging from the start.

---

### Pitfall 12: Hardcoded Region/Zone Without Availability Check

**What goes wrong:**
Spot VM creation fails because chosen zone has no capacity. Script fails without fallback.

**Why it happens:**
Spot capacity varies by zone. A zone that worked yesterday may be unavailable today.

**Prevention:**
1. **Try multiple zones** - Script should retry in different zones
2. **Check availability** - `gcloud compute zones list --filter="region:us-central1"`
3. **For qm-nmr-calc** - Single VM, manual management means manual zone selection is acceptable

**Detection:**
- VM creation fails with resource availability error
- Zone shows as UNAVAILABLE

**Phase to address:** Deployment script phase - handle gracefully or document.

---

### Pitfall 13: Forgetting External IP for Direct Access

**What goes wrong:**
VM created without external IP. Cannot SSH directly or access web interface.

**Why it happens:**
Some examples omit `--address` or use private-only networking for security.

**Prevention:**
```bash
gcloud compute instances create my-vm \
  --network-interface=network-tier=PREMIUM,subnet=default
  # Ephemeral external IP assigned by default
```

Or explicitly:
```bash
gcloud compute addresses create my-static-ip --region=us-central1
gcloud compute instances create my-vm \
  --address=my-static-ip
```

**Detection:**
- `EXTERNAL_IP` column is empty in `gcloud compute instances list`
- Cannot reach VM from internet

**Phase to address:** Deployment script phase - verify external IP allocation.

---

### Pitfall 14: Domain DNS Pointing to Old/Wrong IP

**What goes wrong:**
VM has new ephemeral IP after restart but DNS still points to old IP. HTTPS/Caddy setup fails.

**Why it happens:**
Ephemeral IPs change. DNS TTL may cache old IP.

**Prevention:**
1. **Use static IP** - Reserve IP address, doesn't change across restarts
2. **Update DNS on startup** - Script to update DNS record via Cloud DNS API or external provider
3. **Low TTL during development** - 60-300 seconds instead of 3600+

**Detection:**
- DNS lookup returns wrong IP
- Caddy cannot get Let's Encrypt certificate (HTTP challenge fails)

**Phase to address:** HTTPS/domain setup phase - decide static vs dynamic IP early.

---

## NWChem-Specific Considerations

### Pitfall 15: NWChem Permanent Directory Not on Persistent Disk

**What goes wrong:**
NWChem restart files (`.db`, `.movecs`) written to ephemeral storage. Cannot resume calculations after preemption.

**Why it happens:**
Default NWChem `permanent_dir` may point to local filesystem, not persistent mount.

**Prevention:**
Configure NWChem input:
```
permanent_dir /data/nwchem/permanent
scratch_dir /scratch
```

Mount persistent disk at `/data`, use local SSD (if any) for scratch.

**Detection:**
- Restart attempt fails with "cannot find .db file"
- No `.movecs` files in expected location after restart

**Phase to address:** NWChem integration phase - ensure restart file persistence.

**Sources:**
- [NWChem Restart](https://nwchemgit.github.io/Start_Restart.html)

---

## Phase-Specific Warning Summary

| Phase | Likely Pitfall | Mitigation |
|-------|---------------|------------|
| Infrastructure Setup | Firewall rules not created | Use default VPC or create rules explicitly |
| Infrastructure Setup | Local SSD instead of persistent disk | Always use persistent disk for job data |
| Infrastructure Setup | Wrong OS image | Use COS or Ubuntu (no premium license) |
| Deployment Script | `--preemptible` instead of `--provisioning-model=SPOT` | Use modern Spot VM flags |
| Deployment Script | Termination action DELETE | Always use STOP for stateful workloads |
| Deployment Script | Containers not starting on restart | Use restart policy and startup script |
| HTTPS/Domain Setup | Ephemeral IP changes | Use static IP or update DNS dynamically |
| Operations | No preemption awareness | Add monitoring or accept manual checking |
| NWChem Integration | Restart files on ephemeral storage | Configure permanent_dir on persistent disk |

---

## Cost Gotchas

| Gotcha | Impact | Prevention |
|--------|--------|------------|
| Persistent disk charges continue when VM stopped | ~$0.04/GB/month | Delete data disk if not needed, or accept as backup cost |
| Spot VMs cannot use Committed Use Discounts | No CUD/SUD stacking | Expected - Spot pricing is already 60-91% off |
| Premium OS licensing not discounted | Full OS price even on Spot | Use free OS images |
| Network egress charges | Adds up with large output files | Compress results, use same-region storage |
| Preemption under 1 minute not billed | No refund for partial work | Minor, calculations take hours anyway |

---

## Sources Summary

### Official Documentation (HIGH confidence)
- [GCP Spot VMs](https://cloud.google.com/compute/docs/instances/spot)
- [Create and Use Spot VMs](https://cloud.google.com/compute/docs/instances/create-use-spot)
- [Shutdown Scripts](https://cloud.google.com/compute/docs/shutdownscript)
- [Persistent Disk](https://cloud.google.com/compute/docs/disks/persistent-disks)
- [Firewall Rules](https://cloud.google.com/firewall/docs/firewalls)
- [Container Deployment](https://cloud.google.com/compute/docs/containers/deploying-containers)
- [Artifact Registry Auth](https://cloud.google.com/artifact-registry/docs/docker/authentication)

### Community Resources (MEDIUM confidence)
- [Reliably Executing Shutdown Scripts](https://haggainuchi.com/shutdown.html)
- [Pitfalls in GKE Spot VMs](https://medium.com/google-cloud/pitfalls-to-avoid-when-using-spot-vms-in-gke-for-cost-reduction-c6f42f674c1f)
- [Managing VM Restart in GCP](https://medium.com/google-cloud/managing-vm-restart-in-gcp-all-scenarios-covered-6c861722cf29)

### NWChem Documentation (HIGH confidence)
- [NWChem Start/Restart](https://nwchemgit.github.io/Start_Restart.html)
- [NWChem Architecture](https://nwchemgit.github.io/NWChem-Architecture.html)

---

*Last updated: 2026-02-04*
