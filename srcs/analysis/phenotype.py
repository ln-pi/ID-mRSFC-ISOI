"""Validated optional phenotype/event interfaces; no fabricated observations."""
import pandas as pd

OBSERVATION_COLUMNS=['mouse','date','measure','value']
EVENT_COLUMNS=['mouse','entry_date','event_date','censored']


def read_observations(path, mice):
    frame=pd.read_csv(path)
    if not set(OBSERVATION_COLUMNS)<=set(frame): raise ValueError('Phenotype columns: '+','.join(OBSERVATION_COLUMNS))
    frame['date']=pd.to_datetime(frame['date'],errors='raise')
    frame['value']=pd.to_numeric(frame['value'],errors='raise')
    if frame[OBSERVATION_COLUMNS].isna().any().any() or set(frame.mouse)-set(mice):raise ValueError('Unknown mouse or missing phenotype field')
    if frame.duplicated(['mouse','date','measure']).any():raise ValueError('Duplicate phenotype observation')
    return frame


def align_observations(trajectory, observations, tolerance_days=7):
    """Most recent prior observation per measurement; never use future phenotype."""
    if tolerance_days<0:raise ValueError('tolerance_days must be nonnegative')
    rows=[]
    for _,row in trajectory.iterrows():
        date=pd.to_datetime(row.get('date'),errors='coerce')
        for measure in observations.measure.unique():
            candidates=observations[(observations.mouse==row.mouse)&(observations.measure==measure)]
            candidates=candidates[(candidates.date<=date)&(candidates.date>=date-pd.Timedelta(days=tolerance_days))]
            latest=candidates.sort_values('date').iloc[-1] if len(candidates) else None
            rows.append(dict(row,measure=measure,phenotype_value=latest.value if latest is not None else None,
                phenotype_date=latest.date if latest is not None else None,alignment_status='aligned' if latest is not None else 'no_prior_measurement'))
    return pd.DataFrame(rows)


def read_events(path, mice):
    frame=pd.read_csv(path)
    if not set(EVENT_COLUMNS)<=set(frame):raise ValueError('Event columns: '+','.join(EVENT_COLUMNS))
    for name in ['entry_date','event_date']:frame[name]=pd.to_datetime(frame[name],errors='raise')
    if frame[EVENT_COLUMNS].isna().any().any() or frame.mouse.duplicated().any() or set(frame.mouse)-set(mice):raise ValueError('Invalid event identity or missing field')
    if not frame.censored.isin([0,1,False,True]).all() or (frame.event_date<frame.entry_date).any():raise ValueError('Invalid censoring or event chronology')
    return frame
